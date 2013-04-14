
#ifndef _ABSMTX_H_
#define _ABSMTX_H_
#include <iostream>
#include "vector.h"


/*--------------------------------------------------------------
|
|				Abstract Matrix' Definition
|
--------------------------------------------------------------*/ 
class AbsMtx  
{
public:
	AbsMtx(){};
	virtual ~AbsMtx(){};
public:
	virtual Vec operator*(const Vec &) const = 0;

//Iterative methods for solving linear systems.

	int CG(Vec& x,const Vec& b,double & eps,int & iter,int pcn);
  int GMRERS(Vec& x,const Vec& b,double & eps,int& iter,int pcn,int m);
	int BiCGSTAB(Vec &x, const Vec &b,double &tol,int &max_iter, int pcn);

protected:
	int nrows;
	virtual Vec preconding(const Vec & r,int i=0) const = 0;
};





#endif //_ABSMTX_H_
