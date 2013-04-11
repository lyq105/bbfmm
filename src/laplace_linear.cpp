//==============================================================================
// This program is used to solve laplace equation 
// \Delta u = 0
// with Boundary Element Method.
//==================================================================================*/


#include "mesh.h"
#include "Matrix.h"
#include "elements.h"

/**
 * This is a solver for laplace Equation by using BIEM, linear element.
 *
 * In this class we define the main procedual of 
 * solving this program.
 *
 * @author liyiqiang, 2013
 */

class laplace2d_linear
{
public:
	//! A construct function.
	Laplace2d(){};

	Laplace2d(std::string filename){bmesh.load_mesh(filename);};
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

	struct Surface_element_dof_map
	{
		int number_of_e_dof;
		std::vector<int> e_dof_map;    
	};
	std::vector<Surface_element_dof_map> s_dof_map;

/**
 * Define the boundary condition for this equation.
 */

	int initial_bc();

/**
 * Define the degree of freedom map.
 */

	int init_dofmap();

/**
 * Calculate the Global Matrix and right hand side of this equation.
 */

	int assemble_matrix();

/**
 * Give the main solving procedul of BIEM.
 * we gather all procedual in this function.
 */ 	
	
	int solve();
	
	// Post processing  

/**
 * Get the value at the point which is a field point of region \f$ \Omega \f$
 */
	double u_boundary(int number);
	
/**
 * Get the value at the point which is a field point of region \f$ \Omega \f$
 */
	double u_field(Point& pt);
/**
 *
 */	
	int print_solution(std::string resfname);
	double h0_norm();

private:
	Mesh bmesh;
	Boundary_condition bc;
	FullMtx global_matrix;
	Vec rhs;
	Solution solution;

};

// boundary information
int 
laplace2d_linear::initial_bc()
{
	bc.freedom = bmesh.number_of_surface_elements();
	bc.bctags.resize(bc.freedom);
	bc.bcvalues.resize(bc.freedom);
	//int test = 0 ;
	for (int i = 0; i < bc.freedom; i++){
		Mesh::Surface_element& se = bmesh.find_surface_element(i);
		if(se.se_tag == 1)    // first kind boundary condition
		{
			// 			test++;
			bc.bctags[i] = 1;
			bc.bcvalues[i] = 200;
		}
		if(se.se_tag == 2) // second kind boundary condition
		{ 
			// 			test++;
			bc.bctags[i] = 2;
			bc.bcvalues[i] = 10;
		}
	}
	//	std::cout << "test = " << test << std::endl;
	return 1;
}

int laplace2d_linear::init_dofmap()
{
	for (int i=0;i<s_dof_map.size();++i)
	{
		s_dof_map[i].number_of_e_dof = 2;
		s_dof_map[i].e_dof_map.push_back(bmesh.find_surface_element(i).se_conectivity[0]);
		s_dof_map[i].e_dof_map.push_back(bmesh.find_surface_element(i).se_conectivity[1]);
	}
	return 1;
}

// get the matrix 
int
laplace2d_linear::assemble_matrix()
{
	double Hij = 0;
	double Gij = 0;
	Link2nodeSur l2s;

	int se_number = bmesh.number_of_surface_elements();
	// loop over all surface elements
	for ( int e_id = 0; e_id < se_number; ++e_id )
	{
		// get the element data
		Mesh::Surface_element& se = bmesh.find_surface_element(e_id);
		l2s.initial(bmesh,se);

		Point source[2];
//		source[0] = l2s.to_real(Point(-1,0,0));
//		source[1] = l2s.to_real(Point(1,0,0));

		for (int i = 0; i < 2; ++i){
			for (int j = 0; j < 2; ++j){
				
				int u = s_dof_map[e_id].e_dof_map[i];
				int v = s_dof_map[e_id].e_dof_map[j];

				for (int g_id = 0; g_id < l2s.gauss_int.number_of_gauss_points; ++g_id)
				{
					// TODO: transfer ref point to real plane. 
					Point gauss_pt = l2s.to_real(l2s.gauss_int.p_list[g_id]);
					Hij += l2s.shape(j)*g_kernel_2d(source[i],gauss_pt)
						* l2s.get_det_jacobian(gauss_pt) * l2s.gauss_int.w_list[g_id];
					Gij += l2s.shape(j)*h_kernel_2d(source[i],gauss_pt)
						* l2s.get_det_jacobian(gauss_pt) * l2s.gauss_int.w_list[g_id];
				}	
				
				// assemble the global matrix by using boundary info data.
				if (bc.bctags[u] == 1)  // 1st kind boundary.
				{
					global_matrix[u][v] += -Gij;
					rhs[u] = -Hij*bc.bcvalues[u];
				}
				else if (bc.bctags[u] == 2)
				{
					global_matrix[u][v] += Hij;
					rhs[u] = Gij*bc.bcvalues[u];
				}			
			}
		}
	}
	// Plus I_{\alpha} to global matrix. this part also can be d
	
	return 0;
}

double
laplace2d_linear::solve()
{
	initial_bc();
	init_dofmap();
	assemble_matrix();

	Vec s;
	global_matrix.BiCGSTAB(s,rhs,1e-9,1000,0);

	for (int i = 0; i< bmesh.number_of_s_node; ++i)
	{
		if (bc.bctags[i] == 1)  // 1st kind boundary.
		{
			solution.q[i] = s[i];
			solution.u[i] = bc.bcvalues[i];
		}
		else if (bc.bctags[i] == 2)
		{
			solution.u[i] = s[i];
			solution.q[i] = bc.bcvalues[i];
		}	
	}

}

double g_kernel_2d()
{

}

double h_kernel_2d()
{

}
