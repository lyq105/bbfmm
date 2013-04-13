/**************************************************************
*
*             
*
*This is The linear Algra 's implatement			
*
*
*
*2009.9.16 modified
*
***************************************************************/
#include "vector.h"
#include "absmtx.h"
#include <assert.h>
#include <cmath>
#include <iostream>
using std::cout;
using std::endl;



/*-------------------------------------------------------------
|        
|				 Abstract Matrix 
|
--------------------------------------------------------------*/

int 
AbsMtx::CG(Vec & x,const Vec & b,double & eps,int & iter,int pcn)
{
	assert(nrows == b.size());
	const int maxiter = iter;
	Vec r = b - (*this)*x;
	Vec z = preconding(r,pcn);
	Vec p = z;
	double zr = dot(z,r);
	const double stp = eps /** b.twonorm()*/;

	if (!r.maxnorm())
	{
		eps = 0.0;
		iter = 0;
		return 1;
	}

	for(iter = 0;iter<maxiter;iter++)
	{
		Vec mp =(*this)*p;
		double pap=dot(mp,p);
		double alpha = zr/pap;
		x += alpha*p;
		r -= alpha*mp;
		if (r.twonorm() <= stp)break;
		z = preconding(r,pcn);
		double zrold = zr;
		zr = dot(z,r);
		double beta = zr/zrold;
		p = z+beta*p;
	}

	eps = r.twonorm();
	if (iter == maxiter)return 1;
	else return 0;
}

//return values
//0:
//1-3:	exit normally;
//5: breakdown	
//6: breakdown	


int 
AbsMtx::BiCGSTAB(Vec &x, const Vec &b,double &tol,int &max_iter, int pcn)
{
	
	double resid;
	double rho_1, rho_2, alpha, beta, omega;
	Vec p(b.size()); 	Vec q(b.size());


	double normb = norm(b); 
	Vec r = b - (*this) * x;
	Vec rtilde = r;

	if (normb == 0.0)
		normb = 1;

	if ((resid = norm(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}

	p=r;
	for (int i = 1; i <= max_iter; i++) 
	{
		std::cout << r.twonorm() << std::endl;
		rho_1 = dot(rtilde, r);

		if (rho_1 == 0) {
			tol = norm(r) / normb;
			return 2;
		}
	
		alpha = rho_1 / dot(rtilde, (*this) * p);
		if (fabs(dot(rtilde, (*this) * p))<1e-20)
		{
			tol = norm(r) / normb;
			return 5;
		}
		q = r - alpha * ((*this) * p);
		omega = dot((*this) * q,q) / dot((*this) * q,(*this) * q);
		if (fabs(dot(rtilde, (*this) * q))<1e-20)
		{
			tol = norm(r) / normb;
			return 6;
		}

		if (omega == 0) {
			tol = norm(r) / normb;
			return 3;
		}

		x += alpha * p + omega * q;
		r = q - omega * ((*this) * q);
		
		if ((resid = norm(r) / normb) < tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}
		rho_2=dot(rtilde, r);

		beta = (rho_2/rho_1) * (alpha/omega);
		p = r + beta * (p - omega * ((*this) * p));

	}

	tol = resid;
	return 1;
}



