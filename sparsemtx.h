#ifndef __SPARSE_MATRIX__

#define __SPARSE_MATRIX__



#include "absmtx.h"
#include "fullmtx.h"
#include <iostream>

using std::ostream;

/*-----------------------------------------------------------------
|
|				Sparse Matrix Class'Definition
|
-----------------------------------------------------------------*/
class SparseMtx:public AbsMtx
{
public:
	friend ostream & operator << (ostream &,SparseMtx &);

	//construction and destruction
	SparseMtx(int nrow,int len,double *t,int *c,int *f);
	SparseMtx(int nrow,int len);
	SparseMtx(const SparseMtx &);
	SparseMtx(FullMtx & Fmt);

	~SparseMtx(){delete [] sra;delete [] fnz;delete [] clm;}

public:
	//overload operator
	SparseMtx & operator=(const SparseMtx &);
	Vec operator*(const Vec&) const;
	double operator[](int i)const{return sra[i];}
	double operator()(int i,int j);
	//other
	int getclm(int i) const {return clm[i];}
	int getfnz(int i) const {return fnz[i];}
	double getsra(int i) const {return sra[i];}

	int row()const{return nrows;}
	int getlenth() const {return lenth;}
	int getlenth() {return lenth;}

private:
	int lenth;
	int ncols;
	int * clm;
	int * fnz;
	double * sra;
	Vec preconding(const Vec &,int i = 0) const;

};

#endif //__SPARSE_MATRIX__
