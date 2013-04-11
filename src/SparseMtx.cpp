
#include "SparseMtx.h"
#include <assert.h>
#include <math.h>
/*-------------------------------------------------------
|
|			Sparse Matrix Class
|
--------------------------------------------------------*/

SparseMtx::SparseMtx(int rowsnum,int len,double *t,int *c,int *f)
{
	int i=0;
	nrows = rowsnum;
	lenth = len;
	sra = new double [lenth];
	clm = new int [lenth];
	fnz = new int [nrows+1];

	for ( i=0;i<lenth;i++)
	{
		sra[i]=t[i];
		clm[i]=c[i];
	}
	for( i=0;i<nrows;i++)fnz[i]=f[i];
	fnz[nrows]=lenth;
}

SparseMtx::SparseMtx(int rows,int len)
{
	int i=0;
	nrows = rows;
	lenth = len;
	sra = new double [lenth];
	clm = new int [lenth];
	fnz = new int [nrows+1];

	for ( i=0;i<lenth;i++)
	{
		sra[i]=0;
		clm[i]=0;
	}
	for(  i=0;i<=nrows;i++)fnz[i]=0;
}
SparseMtx::SparseMtx(const SparseMtx &spm)
{
	int i=0;
	nrows = spm.nrows;
	lenth = spm.lenth;
	sra = new double [lenth];
	clm = new int [lenth];
	fnz = new int [nrows+1];

	for ( i=0;i<lenth;i++)
	{
		sra[i]=spm.sra[i];
		clm[i]=spm.clm[i];
	}
	for( i=0;i<=nrows;i++)fnz[i]=spm.fnz[i];
}

SparseMtx::SparseMtx(FullMtx &fmt)
{
	int i,j;
	int nNonzero=0;
	for ( i=0;i<fmt.row();i++)
		for( j=0;j<fmt.row();j++)
			if(fmt[i][j]!=0)nNonzero++;
	nrows=fmt.row();
	lenth=nNonzero;

	sra = new double [lenth];
	clm = new int [lenth];
	fnz = new int [nrows+1];
	int p=0;
	for (i=0;i<fmt.row();i++)
	{
		int k=0;
		for(j=0;j<fmt.row();j++)
			if(fmt[i][j]!=0)
			{
				if(k==0){fnz[i]=p;k=1;}
				sra[p]=fmt[i][j];
				clm[p]=j;
				p++;
			}
	}
	fnz[nrows]=lenth;
}
Vec SparseMtx::operator * (const Vec &v)const
{
	assert(nrows==v.size());
	Vec tm(nrows);
	for (int i=0;i<nrows;i++)
	{
		for (int j=fnz[i];j<fnz[i+1];j++)
		{
			tm[i]+=sra[j]*v[clm[j]];
		}
	}
	return tm;
}
SparseMtx& SparseMtx::operator = (const SparseMtx& spm)
{
	int i=0;
	if(this!=&spm)
	{
		assert(nrows==spm.nrows&&lenth==spm.lenth);
		for ( i=0;i<lenth;i++)
		{
			sra[i]=spm.sra[i];
			clm[i]=spm.clm[i];
		}
		for( i=0;i<=nrows;i++)fnz[i]=spm.fnz[i];
	}
	return *this;
}

double SparseMtx::operator ()(int i,int j)
{
	for(int p=fnz[i];p<fnz[i+1];p++)
		if (clm[p]==j) return sra[p];
	return 0;
}
Vec SparseMtx::preconding(const Vec &r,int precn /* = 0 */)const
{
	int i,j;
	if (precn==0)
	{
		return r;
	}
	else if (precn ==1)
	{
		Vec diag(nrows);
		for( i=0;i<nrows;i++)
		{
			for( j =fnz[i];j<fnz[i+1];j++)
			{
				diag[i]+=sra[j]*sra[j];
			}
			diag[i]=sqrt(diag[i]);
		}
		Vec z(nrows);
		for ( i=0;i<nrows;i++)z[i]=r[i]/diag[i];
		return z;
	}else if(precn==2)
	{
		const double omega =1.2;

		Vec diag(nrows);
		for ( i=0;i<nrows;i++)
		{
			for ( j = fnz[i];j<fnz[i+1];j++)
			{
				if(clm[j]==i)diag[i]+=sra[j]*sra[j];
				else diag[i]+=omega*omega*sra[i]*sra[j];
			}
			diag[i] = sqrt(diag[i]);
		}

		Vec z(nrows);
		for ( i=0;i<nrows;i++)
		{
			double sum = 0;
			for( j =fnz[i];j<fnz[i+1];j++)
				if (clm[j]<i)
				{
					sum+=sra[j]*z[clm[j]];
				}
				z[i]=(r[i]-omega*sum)/diag[i];
		}
		for ( i=nrows-1;i>=0;i--)
		{
			double sum=0;
			for( j=fnz[i];j<fnz[i+1];j++)
			{
				if(clm[j]>i)sum+=sra[j]*z[clm[j]];
			}
			z[i]-=omega*sum/diag[i];
		}
		return z;
	}else{exit(1);}
}
ostream & operator << (ostream & ostm,SparseMtx & spt)
{
	for (int i =0; i<spt.row();i++)
	{
		for (int j =0; j<spt.row();j++)
		{
			ostm<<spt(i,j)<<"\t";
		}
		ostm<<"\n";
	}
	return ostm;
}
