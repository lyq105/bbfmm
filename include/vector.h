#ifndef _VECTOR_H__
#define _VECTOR_H__


#include <iostream>
using std::ostream;

/*--------------------------------------------------------
|
|				Vector Class'Definition
|
---------------------------------------------------------*/
class Vec
{
public:
//construct function

	Vec(int,const double *);
	Vec(int = 1,double =0);
	Vec(const Vec&);

//Destruct function

	~Vec(){delete [] vr;}

//Interface of Vec Class
public:

	Vec &operator=(const Vec&);
	Vec &operator+=(const Vec&);
	Vec &operator-=(const Vec&);
	double &operator [] (int i){return vr[i];}
	double &operator [] (int i)const{return vr[i];}

//
public:

	int size()const {return lenth;}
	double maxnorm() const;
	double twonorm() const;
	bool Setlenth(int len);
//friend function
	friend double dot(const Vec &,const Vec &);
	friend double norm(const Vec&);
	friend double dot(const double*  v1,const double*  v2,int n);

	friend Vec operator+(const Vec&);
	friend Vec operator-(const Vec&);
	friend Vec operator+(const Vec&,const Vec&);
	friend Vec operator-(const Vec&,const Vec&);
	friend Vec operator*(double,const Vec&);
	friend Vec operator*(const Vec&,double);
	friend Vec operator*(const Vec&,const Vec&);
	friend Vec operator/(const Vec&,double);
	friend ostream& operator << (ostream&,const Vec&);

private:
	int lenth;		 //vector's lenth
	double * vr;     //values of vector
};

#endif //_VECTOR_H__
