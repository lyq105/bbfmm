#include "Vector.h"
#include <assert.h>

#include <math.h>
/*-------------------------------------------------------------
|
|				 Vector class
|
--------------------------------------------------------------*/


//construction function
Vec::Vec(int n,const double * abd)
{
	lenth = n;
	vr = new double[lenth];
	for (int i=0;i<lenth;i++){vr[i]=abd[i];}
}
Vec::Vec(const Vec& vf)
{
	lenth=vf.lenth;
	vr = new double[lenth];
	for (int i=0;i<lenth;i++){vr[i]=vf.vr[i];}
}

Vec::Vec(int n,const double valu)
{
	lenth=n;
	vr = new double[lenth];
	for (int i=0;i<lenth;i++){vr[i]=valu;}
}

//operator overload
Vec& Vec::operator =(const Vec& v)
{
	if (this != &v)
	{
		assert(lenth==v.lenth);
		for(int i=0;i<lenth;i++){vr[i]=v.vr[i];}
	}
	return *this;
}
Vec & Vec::operator+=(const Vec& v)
{
	assert(lenth==v.lenth);
	for(int i=0;i<lenth;i++){vr[i]+=v.vr[i];}
	return *this;
}
Vec & Vec::operator-=(const Vec& v)
{
	assert(lenth==v.lenth);
	for(int i=0;i<lenth;i++){vr[i]-=v.vr[i];}
	return *this;
}

//norm
double Vec::maxnorm()const
{
	double nm = fabs(vr[0]);
	for (int i = 1;i<lenth;i++)
	{
		if (nm<=fabs(vr[i]))
		{
			nm=fabs(vr[i]);
		}
	}
	return nm;
}
double Vec::twonorm()const
{
	double norm2 = vr[0]*vr[0];
	for (int i = 1;i<lenth;i++){norm2+=vr[i]*vr[i];}
	return sqrt(norm2);
}

bool Vec::Setlenth(int len)
{
	lenth=len;
	vr = new double[lenth];
	for (int i=0;i<lenth;i++){vr[i]=0;}
	return true;
}



// Vector's Friend Functions

Vec operator+(const Vec & v){return v;}
Vec operator-(const Vec & v){return Vec(v.lenth)-v;}

Vec operator+(const Vec & v1,const Vec & v2)
{
	assert(v1.lenth==v2.lenth);
	Vec sum = v1;
	sum += v2;
	return sum;
}
Vec operator-(const Vec & v1,const Vec & v2)
{
	assert(v1.lenth==v2.lenth);
	Vec ms = v1;
	ms -= v2;
	return ms;
}
Vec operator*(double scalar, const Vec & v)
{
	Vec tm(v.lenth);
	for (int i=0;i<v.lenth;i++)
	{
		tm[i]=scalar*v[i];
	}
	return tm;
}
Vec operator*(const Vec & v,double scalar)
{
	return scalar*v;
}

Vec operator*(const Vec&v1,const Vec&v2)
{
	assert(v1.lenth==v2.lenth);
	Vec tm(v1.lenth);
	for (int i=0;i<v1.lenth;i++)
	{
		tm[i]=v1[i]*v2[i];
	}
	return tm;
}
Vec operator / (const Vec &v,double scalar)
{
	assert(scalar!=0);
	return (1.0/scalar)*v;
}
ostream& operator << (ostream& ostr,const Vec& vt)
{
	for (int i=0;i<vt.lenth;i++)
	{
		ostr<<vt[i]<<"\t";
		if (i%10 == 5){ostr<<"\n";}
	}
	return ostr;
}

double dot(const Vec & v1,const Vec & v2)
{
	assert(v1.lenth == v2.lenth);
	double tm = v1[0]*v2[0];
	for (int i = 1;i<v1.lenth;i++)
	{
		tm+=v1[i]*v2[i];
	}
	return tm;
}
double dot(const double*  v1,const double*  v2,int n)
{
	double tm = 0;
	for (int i = 0;i<n;i++)
	{
		tm+=v1[i]*v2[i];
	}
	return tm;
}
double norm(const Vec& v )
{
	return v.twonorm();
}
