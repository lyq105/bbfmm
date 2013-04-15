#ifndef Laplace2d_H
#define Laplace2d_h

#include "mesh.h"


// Solving a laplace equation by boundary element method.

class Laplace2d
{
public:
	Laplace2d(){};
	Laplace2d(std::string filename){bmesh.load_mesh(filename);bmesh.save(filename + ".nmsh");};
	~Laplace2d(){};

	typedef struct __bc_type
	{
		int freedom;
		std::vector<int> bctags;
		std::vector<double> bcvalues;
	} Boundary_condition;

	typedef struct __solu_type
	{
		std::vector<double> u;
		std::vector<double> q;
	} Solution;

//	int initial_bc(Boundary_condition& bc1);
	int initial_bc();
//	int assemble_matrix();
	int solve();
	int solve_lapack();
	double u_boundary(int number);
	double u_field(double sx0, double sx1);
	int print_solution(std::string resfname);
	double h0_norm();

private:
	Mesh bmesh;
	Boundary_condition bc;
	Solution solution;
};

#endif // Laplace2d_h
