#ifndef __LAPLACE_FEM__
#define __LAPLACE_FEM__


#include <vector>
#include <string>
#include <iomanip>
#include <map>
#include <mkl_lapacke.h>
#include "mesh.h"
#include "elements.h"

class Laplace_fem
{
	public:
		struct Boundary_condition 
		{
			int type;
			int index;
			double value;
		};
		std::vector<Boundary_condition>  bc;

		struct Element_dof_map
		{
			int number_of_e_dof;
			std::vector<int> e_dof_map;    
		};
		std::vector<Element_dof_map> dof_map;

		struct Surface_element_dof_map
		{
			int number_of_e_dof;
			std::vector<int> e_dof_map;    
		};
		std::vector<Surface_element_dof_map> s_dof_map;

		int initial_fem_space(std::string filename);
		int initial_bc();
		int assmble_matrix();
		int apply_bc();
		int solve();
		int plot_result(std::string resfname);
		double rhs_func(Point rpt);
	protected:
		int number_of_dof;
		Mesh mesh;
		double* global_matrix;
		double* rhs;
		//  private:
};


#endif // __LAPLACE_FEM__
