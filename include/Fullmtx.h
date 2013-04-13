
#ifndef __FULL_MATRIX__

#define __FULL_MATRIX__




#include "absmtx.h"
#include <iostream>

/*-----------------------------------------------------------------
|
|				  Dense Matrix 'Definition
|
-----------------------------------------------------------------*/
class FullMtx:public AbsMtx
{
public:
//construct
	FullMtx(int n,int m,double**);
	FullMtx(int n=1,int m=1,double t = 0);
	FullMtx(const FullMtx &);

//destruct
	~FullMtx()
	{
		for (int i=0;i<nrows;i++)
		{
			delete[] mx[i];
		}
		delete mx;
	}
//interface
public:
	friend ostream & operator << (ostream &,FullMtx &);

// overloaded operators
	FullMtx & operator=(const FullMtx &);
	FullMtx & operator+=(const FullMtx &);
	FullMtx & operator-=(const FullMtx &);
	FullMtx  operator +(const FullMtx &);
	FullMtx  operator -(const FullMtx &);


	FullMtx & operator+();
	FullMtx  operator-();
	inline double* operator[](int i){return mx[i];}
	inline double* operator[](int i)const{return mx[i];}
	Vec operator*(const Vec&) const;
	Vec operator*(const Vec&);
// Gauss Elimination Method
	int GAUSSELIMINATE(Vec& x,Vec& b) const;
// row and column
	inline int row(){return nrows;}
	inline int column(){return ncols;}
	int insert(int i,int j,double val);
	int Setsize(int n,int m);
private:
	int ncols;
	double ** mx;
	Vec preconding(const Vec & r,int i=0) const;
};

#endif  //__FULL_MATRIX__
