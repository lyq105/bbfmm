#ifndef _SPARSEMTXEX_H_
#define _SPARSEMTXEX_H_

#include "vector.h"
#include "matrix.h"
#include <vector>
using std::vector;


//稀疏方阵
class SparseMtxEx:public AbsMtx 
{
public:
	friend ostream & operator << (ostream &,SparseMtxEx &);
	//construction and destruction 
	SparseMtxEx(int nrow,int len,double *t,int *c,int *f);
	SparseMtxEx(int nrow = 1);
	SparseMtxEx(const SparseMtxEx &);
	SparseMtxEx(SparseMtx &);
	SparseMtxEx(FullMtx & Fmt);
	//
	~SparseMtxEx(){delete [] fnz;}
public:
	//overload operator
	SparseMtxEx & operator=(const SparseMtxEx &);
	Vec operator*(const Vec&) const;
	double operator[](int i)const{return sra[i];}
	double operator()(int i,int j);

	bool Setrow(int nrow);

	//other
	int getclm(int i) const {return clm[i];}
	int getfnz(int i) const {return fnz[i];}
	double getsra(int i) const {return sra[i];}

	int row()const{return nrows;}
	int getlenth(){return lenth;}
	int getlenth() const {return lenth;}

	int insert(int i,int j,double value);

private:
	int lenth;			//非零元素个数
	int* fnz;			//每行第一个非零元素在sra中的位置
	vector<int> clm;    //sra中相应位置的元素的列号
	vector<double> sra; //非零元素的值
	Vec preconding(const Vec &,int i = 0) const;
};


#endif //_SPARSEMTXEX_H_
