#include "FullMtx.h"
#include <assert.h>
#include <math.h>

/*-------------------------------------------------------------
|
|				Full matrix
|
--------------------------------------------------------------*/

FullMtx::FullMtx(int n,int m,double** dbp)
{
	nrows = n;
	ncols = m;
	mx = new double*[nrows];
	for (int i=0;i<nrows;i++)
	{
		mx[i]=new double[ncols];
		for(int j=0;j<ncols;j++)mx[i][j]=dbp[i][j];
	}
}

FullMtx::FullMtx(int n,int m,double t /* = 0 */)
{
	nrows = n;
	ncols = m;
	mx = new double*[nrows];
	for (int i=0;i<nrows;i++)
	{
		mx[i]=new double[ncols];
		for(int j=0;j<ncols;j++)mx[i][j]=t;
	}
}
FullMtx::FullMtx(const FullMtx &fmt)
{
	nrows = fmt.nrows;
	ncols = fmt.ncols;
	mx = new double*[nrows];
	for (int i=0;i<nrows;i++)
	{
		mx[i] = new double [ncols];
		for (int j = 0;j<ncols;j++){mx[i][j]=fmt.mx[i][j];}
	}
}

FullMtx& FullMtx::operator = (const FullMtx& fmt)
{
	if (this != & fmt)
	{
		assert(nrows==fmt.nrows&&ncols==fmt.ncols);
		for (int i = 0;i < nrows;i++)
			for(int j = 0;j < ncols; j++)
				mx[i][j]=fmt.mx[i][j];
	}
	return *this;
}

FullMtx& FullMtx::operator += (const FullMtx& fmt)
{
	assert(nrows==fmt.nrows&&ncols==fmt.ncols);
	for (int i = 0;i < nrows;i++)
		for(int j = 0;j < ncols; j++)
			mx[i][j]+=fmt.mx[i][j];
	return *this;
}

FullMtx& FullMtx::operator -=(const FullMtx& fmt)
{
	assert(nrows==fmt.nrows&&ncols==fmt.ncols);
	for (int i = 0;i < nrows;i++)
		for(int j = 0;j < ncols; j++)
			mx[i][j]-=fmt.mx[i][j];
	return *this;
}
FullMtx& FullMtx::operator +()
{
	return *this;
}

FullMtx FullMtx::operator -()
{
	FullMtx zero(nrows,ncols);
	return zero - *this;
}

FullMtx FullMtx::operator +(const FullMtx & fmt)
{
	FullMtx sum =*this;
	sum +=fmt;
	return sum;
}
FullMtx FullMtx::operator -(const FullMtx & fmt)
{
	FullMtx sum =*this;
	sum -=fmt;
	return sum;
}
Vec FullMtx::operator*(const Vec& vt)
{
	assert(ncols==vt.size());
	Vec tm(nrows);
	for (int i=0;i<nrows;i++)
		for(int j=0;j<ncols;j++)
			tm[i]+=mx[i][j]*vt[j];
	return tm;
}
Vec FullMtx::operator*(const Vec& vt) const
{
	assert(ncols==vt.size());
	Vec tm(nrows);
	for (int i=0;i<nrows;i++)
		for(int j=0;j<ncols;j++)
			tm[i]+=mx[i][j]*vt[j];
	return tm;
}


Vec FullMtx::preconding(const Vec& r,int precn)const
{
	int i=0;
	if (precn==0)
	{
		return r;
	}
	else if(precn ==1)
	{
		Vec z(nrows);
		for ( i = 0; i<nrows;i++)
		{
			z[i]=r[i]/mx[i][i];
		}
		return z;
	}
	else if (precn ==2)
	{
		const double omega=1.2;
		Vec z(nrows);
		for ( i = 0;i<nrows;i++)
		{
			double sum = 0;
			for(int j = 0;j<i;j++)sum+=mx[i][j]*z[j];
			z[i] = (r[i]-omega*sum)/mx[i][i];
		}
		for( i = nrows-1;i>=0;i--)
		{
			double sum=0;
			for (int j = i+1;j<nrows;j++)sum+=mx[i][j]*z[j];
			z[i]-= omega*sum/mx[i][i];
		}
		return z;
	}else {exit(1);}
}


// Friend Methods
ostream & operator << (ostream & ostm,FullMtx & fmt)
{
	for (int i =0; i<fmt.row();i++)
	{
		for (int j =0; j<fmt.column();j++)
		{
			ostm<<fmt[i][j]<<"\t";
		}
		ostm<<"\n";
	}
	return ostm;
}


int FullMtx::insert(int i,int j,double val)
{
	assert(i<nrows&&j<ncols);
	mx[i][j]=val;
	return 1;
}
int FullMtx::Setsize(int n,int m)
{
	nrows = n;
	ncols = m;
	mx = new double*[nrows];
	for (int i=0;i<nrows;i++)
	{
		mx[i]=new double[ncols];
		for(int j=0;j<ncols;j++)mx[i][j]=0;
	}
	return 0;
}

int 
FullMtx::GAUSSELIMINATE(Vec& x,Vec& b) const
{
	int i;
	int k;
	for(k=0;k<ncols-1;k++)
	{
		// 选主元
		double bmax=0.0;
		int ik;
		for(i=k;i<ncols;i++)
		{
			if(bmax<fabs(mx[i][k]))
			{
				bmax=fabs(mx[i][k]);
				ik=i;
			}
		}
		if(bmax<1.0e-10) 
		{
			return 1;
		}
		// 交换第ik行和第k行的元素
		if(ik!=k)
		{
			double t;
			for(i=k;i<ncols;i++)
			{
				t=mx[ik][i];
				mx[ik][i]=mx[k][i];
				mx[k][i]=t;
			}
			t=b[ik];
			b[ik]=b[k];
			b[k]=t;
		}
		// 消元
		for(i=k+1;i<ncols;i++)
		{
			if(mx[i][k]!=0)
			{
				double lk=mx[i][k]/mx[k][k];
				int j;
				for(j=k+1;j<ncols;j++)
				{
					mx[i][j]=mx[i][j]-lk*mx[k][j];
				}
				b[i]=b[i]-lk*b[k];
			}
		}
	}
	if(fabs(mx[ncols-1][ncols-1])<1.0e-10)
	{
		return 1;
	}

	// 消元法结束后开始回代
	x[ncols-1]=b[ncols-1]/mx[ncols-1][ncols-1];
	for(i=ncols-2;i>=0;i--)
	{
		double s=0.0;
		int j;
		for(j=i+1;j<ncols;j++)
		{
			s=s+mx[i][j]*b[j];
		}
		x[i]=(b[i]-s)/mx[i][i];
	}
	return 0;
};
