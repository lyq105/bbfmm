/// using constant boundary element method to solve coupled radiative
//  heat transfer problem by Newton's method.

#ifndef __RADIATIVE_H
#define __RADIATIVE_H

#include "lapace2d.h"
#include <map>
#include <vector>
#include <mkl_lapacke.h>

class CoupledRadiative2d: public Laplace2d
{
public:
	CoupledRadiative2d(){};
	CoupledRadiative2d(std::string filename);
	~CoupledRadiative2d(){};

	// Define the radiative boundary
	typedef struct __rbc_type
	{
		std::vector<int> rbc_tag;
		int rbc_freedom;
		std::vector<int> mapping;
		std::vector<double> u;
		std::vector<double> q;
	} Radiative_boundary;

	struct Element_dof_map{
		int number_of_dof;
		std::vector<int> totle_dof; 
	};

	int initial_bc();
	int load_mesh(std::string filename);
	int solve_newton();
	int cal_rhs();
	int cal_jacobian_matrix();
	int find_number_of_freedom(){return number_of_freedom;}
private:
	Mesh bmesh;
	int number_of_freedom;
	double * jacobian_matrix;
	double * rhs;
	Boundary_condition bc;
	Radiative_boundary rbc;
	Solution solution;
	std::vector<Element_dof_map> dof_map;
};

CoupledRadiative2d::CoupledRadiative2d(std::string filename)
{
	bmesh.load_mesh(filename);
}

int CoupledRadiative2d::load_mesh(std::string filename)
{
	return  bmesh.load_mesh(filename);
};

int CoupledRadiative2d::initial_bc()
{
	bc.freedom = bmesh.number_of_surface_elements();
	bc.bctags.resize(bc.freedom);
	bc.bcvalues.resize(bc.freedom);

	// Degree of freedom mapping's length is number of surface elements.

	dof_map.resize(bmesh.number_of_surface_elements());

	int number_of_radiative_node = 0 ;
	for (int i = 0; i < bc.freedom; i++)
	{	
		Mesh::Surface_element& se = bmesh.find_surface_element(i);
		if(se.se_tag == 1)    // first kind boundary condition
		{
			bc.bctags[i] = 1;
			bc.bcvalues[i] = 100;
			dof_map[i].number_of_dof = 1;
			dof_map[i].totle_dof.push_back(i);
		}
		if(se.se_tag == 2) // second kind boundary condition
		{ 
			bc.bctags[i] = 2;
			bc.bcvalues[i] = 200;
			dof_map[i].number_of_dof = 1;
			dof_map[i].totle_dof.push_back(i);
		}
		if (se.se_tag == 3)
		{
			rbc.rbc_tag[i] = se.se_mat;
			dof_map[i].number_of_dof = 2;
			dof_map[i].totle_dof.push_back(i);
			dof_map[i].totle_dof.push_back(number_of_radiative_node + bmesh.number_of_surface_elements());
			number_of_radiative_node ++;
		}
	}
	rbc.rbc_freedom = number_of_radiative_node;
	number_of_freedom = bmesh.number_of_surface_elements() + number_of_radiative_node;

	return 1;
}


int CoupledRadiative2d::cal_rhs()
{
	const double sigma = 5.670373E-8 ;
	const double epsilong = 0.1;

	int i,j;
	double hij,gij,x,y;

	const double pi = 3.14159265;
	double gaussnum = 4;
	double gausspt[4] = { 0.86113631, -0.86113631, 0.33998014, -0.33998014 };
	double gausswt[4] = {0.34785485,0.34785485,0.65214515,0.65214515};

	//double gaussnum = 7;
	//double gausspt[7] = {0, 0.94910791, -0.94910791, 0.74153119, -0.74153119, 0.40584515,-0.40584515 };
	//double gausswt[7] = {0.41795918,0.12948497,0.12948497,0.27970539,0.27970539,0.38183005,0.38183005};

	int freedom = number_of_freedom;


	for (i = 0; i < freedom; i++){rhs[i] = 0;}


	for(i = 0; i< freedom; i++) 
	{
		Mesh::Surface_element& source_e = bmesh.find_surface_element(i);
		Mesh::Point& s1 = bmesh.find_point(source_e.se_conectivity[0]);
		Mesh::Point& s2 = bmesh.find_point(source_e.se_conectivity[1]);
		// find out the source point
		double source_x0 = 0.5 * (s1.x[0] + s2.x[0]);
		double source_x1 = 0.5 * (s1.x[1] + s2.x[1]);

		double unxi =  source_e.se_outer_norm[0];
		double unyi =  source_e.se_outer_norm[1];

		
		if (bc.bctags[i] != 3)  // Heat conduction equation.
		{
			for (j = 0; j < freedom; j++) 
			{
				// find out the number of 2 ends of a cell.
				Mesh::Surface_element& se = bmesh.find_surface_element(j);
				Mesh::Point& p1 = bmesh.find_point(se.se_conectivity[0]);
				Mesh::Point& p2 = bmesh.find_point(se.se_conectivity[1]);

				// calculate the length of a cell.
				double gama_j = sqrt((p1.x[0] - p2.x[0])*(p1.x[0] - p2.x[0]) 
					+ (p1.x[1] - p2.x[1])*(p1.x[1] - p2.x[1]));
				// calculate the outer unit normal vector of the cell.
				double unx =  se.se_outer_norm[0];
				double uny =  se.se_outer_norm[1];
				if(i == j) // Diagonal element of matrix which have singular integration.
				{
					hij = 1.0/2;
					gij = gama_j / (2*pi) * (log(2.0/gama_j) + 1);
				}
				else
				{
					// Use gauss integral to calculate
					// \f[
					//  {h_{ij}} = \int_{{\Gamma _e}}^{} {F(x,y)} d{\Gamma _e} 
					//     = \sum\limits_{}^{} {F(x,{y_i})*0.5*d*{w_i}}
					//  {g_{ij}} = \int_{{\Gamma _e}}^{} {G(x,y)} d{\Gamma _e}
					//    = \sum\limits_{}^{} {G(x,{y_i})*0.5*d*{w_i}}
					// \f]
					hij = 0;
					gij = 0;
					for(int gindex = 0; gindex< gaussnum; gindex++)
					{
						x = 0.5*(p1.x[0] - p2.x[0])*gausspt[gindex] + 0.5*(p1.x[0] + p2.x[0]);
						y = 0.5*(p1.x[1] - p2.x[1])*gausspt[gindex] + 0.5*(p1.x[1] + p2.x[1]);
						double rr = ((x-source_x0)*(x-source_x0) + (y-source_x1)*(y-source_x1));
						hij += 1.0/ rr * gausswt[gindex];
						gij += - 0.5*log((rr)) *gausswt[gindex];
					}
					double d = (x-source_x0)*unx + (y-source_x1)*uny;
					hij *= -1.0/(4*pi)*gama_j*d;
					gij *= 1.0/(4*pi)*gama_j;
				}
				// Assemble the right hand side vector;
				
				rhs[i] += -hij*solution.u[j] + gij*solution.q[j];
			}
		}
		else   // Surface radiative heat boundary equation.
		{
			for (j = 0; j < freedom; j++) 
			{
				if (rbc.rbc_tag[i] == rbc.rbc_tag[j])   // for elements in the same circle.
				{
					// find out the number of 2 ends of a cell.
					Mesh::Surface_element& se = bmesh.find_surface_element(j);
					Mesh::Point& p1 = bmesh.find_point(se.se_conectivity[0]);
					Mesh::Point& p2 = bmesh.find_point(se.se_conectivity[1]);

					// calculate the length of a cell.
					double gama_j = sqrt((p1.x[0] - p2.x[0])*(p1.x[0] - p2.x[0]) 
						+ (p1.x[1] - p2.x[1])*(p1.x[1] - p2.x[1]));
					// calculate the outer unit normal vector of the cell.
					double unx =  se.se_outer_norm[0];
					double uny =  se.se_outer_norm[1];

					if(i == j) // Diagonal element of matrix which have singular integration.
					{
						hij = epsilong*sigma*pow(solution.u[j],4);
						gij = 1;
					}
					else
					{
						hij = 0;
						gij = 0;
						double mtemp = 0;
						for(int gindex = 0; gindex< gaussnum; gindex++)
						{
							double yy0 = 0.5*(p1.x[0] - p2.x[0])*gausspt[gindex] + 0.5*(p1.x[0] + p2.x[0]);
							double yy1 = 0.5*(p1.x[1] - p2.x[1])*gausspt[gindex] + 0.5*(p1.x[1] + p2.x[1]);
							double fz = (unx*(source_x0 - yy0)+unxi*(source_x1 - yy1))*(uny*(yy0-source_x0)+unyi*(yy1-source_x1));
							double rr = ((yy0-source_x0)*(yy0-source_x0) + (yy1-source_x1)*(yy1-source_x1));
							mtemp += fabs(fz)/(2*rr*sqrt(rr))*gausswt[gindex];
						}

						hij = -epsilong*(1/epsilong -1)*gama_j*mtemp*0.5;
						gij = - epsilong*sigma*pow(solution.u[j],4)*(gama_j*0.5*mtemp);
					}
					// Assemble the matrix and right hand side vector;

					rhs[i] += -hij*solution.u[j] + gij*solution.q[j];

				}
			}
		}
	} // Loop for source point .

	return 0;

}
int CoupledRadiative2d::cal_jacobian_matrix()
{
	const double sigma = 5.670373E-8 ;
	const double epsilong = 0.1;

	int i,j;
	double hij,gij,x,y;

	const double pi = 3.14159265;
	double gaussnum = 4;
	double gausspt[4] = { 0.86113631, -0.86113631, 0.33998014, -0.33998014 };
	double gausswt[4] = {0.34785485,0.34785485,0.65214515,0.65214515};

	//double gaussnum = 7;
	//double gausspt[7] = {0, 0.94910791, -0.94910791, 0.74153119, -0.74153119, 0.40584515,-0.40584515 };
	//double gausswt[7] = {0.41795918,0.12948497,0.12948497,0.27970539,0.27970539,0.38183005,0.38183005};

	int freedom = number_of_freedom;

	
	for (i = 0; i < freedom * freedom; i++){jacobian_matrix[i] = 0;}


	for(i = 0; i< freedom; i++) 
	{
		Mesh::Surface_element& source_e = bmesh.find_surface_element(i);
		Mesh::Point& s1 = bmesh.find_point(source_e.se_conectivity[0]);
		Mesh::Point& s2 = bmesh.find_point(source_e.se_conectivity[1]);
		// find out the source point
		double source_x0 = 0.5 * (s1.x[0] + s2.x[0]);
		double source_x1 = 0.5 * (s1.x[1] + s2.x[1]);

		double unxi =  source_e.se_outer_norm[0];
		double unyi =  source_e.se_outer_norm[1];

		// Heat conduction equation.
		if (bc.bctags[i] != 3)
		{
			for (j = 0; j < freedom; j++) 
			{
				// find out the number of 2 ends of a cell.
				Mesh::Surface_element& se = bmesh.find_surface_element(j);
				Mesh::Point& p1 = bmesh.find_point(se.se_conectivity[0]);
				Mesh::Point& p2 = bmesh.find_point(se.se_conectivity[1]);

				// calculate the length of a cell.
				double gama_j = sqrt((p1.x[0] - p2.x[0])*(p1.x[0] - p2.x[0]) 
					+ (p1.x[1] - p2.x[1])*(p1.x[1] - p2.x[1]));
				// calculate the outer unit normal vector of the cell.
				double unx =  se.se_outer_norm[0];
				double uny =  se.se_outer_norm[1];
				if(i == j) // Diagonal element of matrix which have singular integration.
				{
					hij = 1.0/2;
					gij = gama_j / (2*pi) * (log(2.0/gama_j) + 1);
				}
				else
				{
					// Use gauss integral to calculate
					// \f[
					//  {h_{ij}} = \int_{{\Gamma _e}}^{} {F(x,y)} d{\Gamma _e} 
					//     = \sum\limits_{}^{} {F(x,{y_i})*0.5*d*{w_i}}
					//  {g_{ij}} = \int_{{\Gamma _e}}^{} {G(x,y)} d{\Gamma _e}
					//    = \sum\limits_{}^{} {G(x,{y_i})*0.5*d*{w_i}}
					// \f]
					hij = 0;
					gij = 0;
					for(int gindex = 0; gindex< gaussnum; gindex++)
					{
						x = 0.5*(p1.x[0] - p2.x[0])*gausspt[gindex] + 0.5*(p1.x[0] + p2.x[0]);
						y = 0.5*(p1.x[1] - p2.x[1])*gausspt[gindex] + 0.5*(p1.x[1] + p2.x[1]);
						double rr = ((x-source_x0)*(x-source_x0) + (y-source_x1)*(y-source_x1));
						hij += 1.0/ rr * gausswt[gindex];
						gij += - 0.5*log((rr)) *gausswt[gindex];
					}
					double d = (x-source_x0)*unx + (y-source_x1)*uny;
					hij *= -1.0/(4*pi)*gama_j*d;
					gij *= 1.0/(4*pi)*gama_j;
				}
				// Assemble the matrix and right hand side vector;
				if(bc.bctags[j] == 1)
				{
					int u = dof_map[i].totle_dof[0];
					int v = dof_map[j].totle_dof[0];
					jacobian_matrix[u + v *freedom] = -gij;
				}
				else if (bc.bctags[j] == 2)
				{
					int u = dof_map[i].totle_dof[0];
					int v = dof_map[j].totle_dof[0];
					jacobian_matrix[u + v *freedom] = hij;
				}
				else if (bc.bctags[j] == 3)
				{
					int u = dof_map[i].totle_dof[0];
					int v = dof_map[j].totle_dof[0];
					int w = dof_map[j].totle_dof[1];
					jacobian_matrix[u + v *freedom] = hij;
					jacobian_matrix[u + w *freedom] = -gij;
				}
			}
		}
		else   // Surface radiative heat boundary equation.
		{
			for (j = 0; j < freedom; j++) 
			{
				if (rbc.rbc_tag[i] == rbc.rbc_tag[j])   // for elements in the same circle.
				{
					// find out the number of 2 ends of a cell.
					Mesh::Surface_element& se = bmesh.find_surface_element(j);
					Mesh::Point& p1 = bmesh.find_point(se.se_conectivity[0]);
					Mesh::Point& p2 = bmesh.find_point(se.se_conectivity[1]);

					// calculate the length of a cell.
					double gama_j = sqrt((p1.x[0] - p2.x[0])*(p1.x[0] - p2.x[0]) 
						+ (p1.x[1] - p2.x[1])*(p1.x[1] - p2.x[1]));
					// calculate the outer unit normal vector of the cell.
					double unx =  se.se_outer_norm[0];
					double uny =  se.se_outer_norm[1];

					if(i == j) // Diagonal element of matrix which have singular integration.
					{
						hij = 4*epsilong*sigma*pow(solution.u[j],3);
						gij = 1;
					}
					else
					{
						hij = 0;
						gij = 0;
						double mtemp = 0;
						for(int gindex = 0; gindex< gaussnum; gindex++)
						{
							double yy0 = 0.5*(p1.x[0] - p2.x[0])*gausspt[gindex] + 0.5*(p1.x[0] + p2.x[0]);
							double yy1 = 0.5*(p1.x[1] - p2.x[1])*gausspt[gindex] + 0.5*(p1.x[1] + p2.x[1]);
							double fz = (unx*(source_x0 - yy0)+unxi*(source_x1 - yy1))*(uny*(yy0-source_x0)+unyi*(yy1-source_x1));
							double rr = ((yy0-source_x0)*(yy0-source_x0) + (yy1-source_x1)*(yy1-source_x1));
							mtemp += fabs(fz)/(2*rr*sqrt(rr))*gausswt[gindex];
						}

						hij = -epsilong*(1/epsilong -1)*gama_j*mtemp*0.5;
						gij = - 4*epsilong*sigma*pow(solution.u[j],3)*(gama_j*0.5*mtemp);
					}
					// Assemble the matrix and right hand side vector;
					int u = dof_map[i].totle_dof[0];
					int v = dof_map[j].totle_dof[0];

					jacobian_matrix[i + j *freedom] += hij;
					jacobian_matrix[i + j *freedom + bmesh.number_of_surface_elements()] += -gij;

				}
			}
		}
	} // Loop for source point .

	return 0;
}

int CoupledRadiative2d::solve_newton()
{
	// using newton's method to solve this problem.
	// initial 

	int freedom = number_of_freedom;

	jacobian_matrix = new double [number_of_freedom*number_of_freedom];
	rhs = new double [number_of_freedom*number_of_freedom];

	// initial guess of solution.

	for (int i = 0; i < freedom; i++) 
	{
		//		std::cout << rhs[i]<<std::endl;
		if (bc.bctags[i] == 1)
		{
			solution.u[i] = bc.bcvalues[i];
			solution.q[i] = 0;
		}
		if (bc.bctags[i] == 2)
		{
			solution.u[i] = 0;
			solution.q[i] = bc.bcvalues[i];
		}
		if (bc.bctags[i] == 3)
		{
			solution.u[i] = 0;
			solution.q[i] = 0;
		}
	}
	// iterative procedure
	double newton_residue = 0;
	do{
		// calculate the equation systems.
		cal_jacobian_matrix();
		cal_rhs();

		// solve equation.
		int* ipiv = new int[freedom];
		int n = freedom;
		int nrhs = 1;
		int lda = freedom;
		int ldb = freedom;

		int info = LAPACKE_dgesv(LAPACK_COL_MAJOR,n,nrhs,jacobian_matrix,lda,ipiv,rhs,ldb);
		delete [] ipiv;

		if (info != 0)
		{
			std::cout << "Calling lapack_dgesv failed!!"  << std::endl; exit(100);
		}

		// update the solution. 

		for (int i = 0; i < freedom; i++) 
		{
			//		std::cout << rhs[i]<<std::endl;
			if (bc.bctags[i] == 1)
			{
				solution.q[i] = rhs[i];
			}
			if (bc.bctags[i] == 2)
			{
				solution.u[i] = rhs[i];
			}
			if (bc.bctags[i] == 3)
			{
				solution.u[i] = rhs[i];
				solution.q[i] = rhs[i + bmesh.number_of_surface_elements()];
			}
		}

		// calculate the residue of Newton's iteration.
		newton_residue = 0;
		for (int i = 0; i < freedom; ++i)
		{
			newton_residue += rhs[i]*rhs[i];
		}
		newton_residue = sqrt(newton_residue);
	} while(newton_residue > 1e-9);// end of newton's iteration. 

	//  free memory
	delete [] jacobian_matrix;
	delete [] rhs;
};

double g_kernel_2d(double x1, double x2,double y1,double y2)
{
	return 1;
}
double f_kernel_2d(double x1, double x2,double y1,double y2)
{
	return 1;
}

double g_kernel_3d(double x1, double x2,double y1,double y2)
{
	return 1;
}
double f_kernel_3d(double x1, double x2,double y1,double y2)
{
	return 1;
}

typedef Mesh::Point Point;
double r_kernel_2d(Point x, Point y, Point nx, Point ny)
{
	return 1;
}
double r_kernel_3d(Point x, Point y, Point nx, Point ny)
{
	return 1;
}
#endif
