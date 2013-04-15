#include <vector>
#include <string>
#include <iomanip>
#include <map>
#include <mkl_lapacke.h>
#include "mesh.h"
#include "elements.h"
#include "laplace_fem.h"

int Laplace_fem::initial_fem_space(std::string filename)
{
	// load mesh;
	mesh.load_mesh(filename);
	// initial degree of freedom map;
	for (int i = 0; i < mesh.number_of_volume_elements(); i++) 
	{
		Laplace_fem::Element_dof_map edof;
		Mesh::Volume_element& ve = mesh.find_volume_element(i);
		edof.number_of_e_dof = 3;

		for (int j=0; j< edof.number_of_e_dof; ++j)
		{
			edof.e_dof_map.push_back(ve.ve_conectivity[j]);
		}    
		dof_map.push_back(edof);
	}

	for (int i = 0; i < mesh.number_of_surface_elements(); i++) 
	{
		Laplace_fem::Surface_element_dof_map sedof;
		Mesh::Surface_element& se = mesh.find_surface_element(i);
		sedof.number_of_e_dof = 2;

		for (int j=0; j< sedof.number_of_e_dof; ++j)
		{
			sedof.e_dof_map.push_back(se.se_conectivity[j]);
		}    
		s_dof_map.push_back(sedof);
	}
	number_of_dof ==  mesh.number_of_point();
	return 0;
}

int Laplace_fem::assmble_matrix()
{
	// loop of all elements
	for ( int e_id = 0; e_id < number_of_dof; ++ e_id )
	{
		// calculate the elements matrix.
		// get the element data
		Mesh::Volume_element& ve = mesh.find_volume_element(e_id);

		Tri3nodeVol t3v;
		t3v.initial(mesh,ve);
		double ke[3][3] = {0};
		double fe[3] = {0};
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				for (int g_id = 0; g_id < t3v.gauss_int.number_of_gauss_points; ++g_id)
				{
					ke[i][j] += ( t3v.dshape_real(i,0,t3v.gauss_int.p_list[g_id])*t3v.dshape_real(j,0,t3v.gauss_int.p_list[g_id])
							+ t3v.dshape_real(i,1,t3v.gauss_int.p_list[g_id])*t3v.dshape_real(j,1,t3v.gauss_int.p_list[g_id]) )
						* t3v.gauss_int.w_list[g_id]*t3v.get_det_jacobian();
				}
				int u = dof_map[e_id].e_dof_map[i];
				int v = dof_map[e_id].e_dof_map[j];
				global_matrix[u + v * number_of_dof] = ke[i][j];
			}
		}

		for (int i = 0; i < 3; ++i)
		{
			for (int g_id = 0; g_id < t3v.gauss_int.number_of_gauss_points; ++g_id)
			{
				fe[i] += t3v.shape(i,t3v.gauss_int.p_list[g_id])* rhs_func(t3v.gauss_int.p_list[g_id])
					* t3v.gauss_int.w_list[g_id]*t3v.get_det_jacobian();
			}
			int u = dof_map[e_id].e_dof_map[i];
			rhs[u] = fe[i];
		}
		// assemble element matrix to global matrix.
	}
	return 0;
}

int Laplace_fem::apply_bc()
{
	for (int b_id =0; b_id<bc.size();++b_id)
	{
		if (bc[b_id].type == 1) // first kind bc;
		{
			int index = bc[b_id].index;
			for (int f_id = 0; f_id < number_of_dof; ++f_id )
			{
				rhs[f_id] -= global_matrix[f_id + number_of_dof* index] * bc[b_id].value;
				global_matrix[index + number_of_dof* f_id] = 0.;
				global_matrix[f_id + number_of_dof* index] = 0.;
			}
			global_matrix[index + number_of_dof*index] = 1.0;
			rhs[index] = bc[b_id].value;
		}
		if (bc[b_id].type == 2) // second kind bc;
		{
			const Mesh::Surface_element& se = mesh.find_surface_element(bc[b_id].index);
			Link2nodeSur l2s;
			l2s.initial(mesh,se);
			double fe[2] = {0};
			for (int i = 0; i<2; ++i)
			{
				for (int g_id = 0; g_id<l2s.gauss_int.number_of_gauss_points; ++ g_id)
				{
					Point& real_pt = l2s.to_real(l2s.gauss_int.p_list[g_id]);
					fe[i] += l2s.get_det_jacobian()*l2s.gauss_int.w_list[g_id] 
					*real_pt.x[i];// given flux on the boundary.
				}
				int u = s_dof_map[bc[b_id].index].e_dof_map[i];
				rhs[u] += fe[i];
			}
		}
	}
	return 1;
}

int Laplace_fem::solve()
{
	
	initial_bc();
	assmble_matrix();
	apply_bc();
	
	int freedom = number_of_dof;
	int* ipiv = new int[freedom];
	int n = freedom;
	int nrhs = 1;
	int lda = freedom;
	int ldb = freedom;

	int info = LAPACKE_dgesv(LAPACK_COL_MAJOR,n,nrhs,global_matrix,lda,ipiv,rhs,ldb);
	delete [] ipiv;

	if (info == 0){
		std::cout << " Successfully call Lapack_dgesv!!"  << std::endl;
	}
	else{
		std::cout << "Calling lapack_dgesv failed!!"  << std::endl;
		exit(100);
	}
//	plot_result();
	return 1;
}

int Laplace_fem::plot_result(std::string resfname)
{
	std::ofstream ofile(resfname.c_str(),std::ios::out);

	ofile << "zone f=fepoint,e=" 
		<< mesh.number_of_volume_elements()
		<< ", n= " 
		<< mesh.number_of_point()
		<< ",et=triangle"
		<< std::endl;


	for (int i = 0; i< mesh.number_of_point(); i++)
	{
		Mesh::Point& p0 = mesh.find_point(i);

		ofile << std::setw(18) << std::setprecision(8)
			<< p0.x[0] << " "
			<< p0.x[1] << " "
			<< rhs[i]
				<< std::endl;
	}
	for (int i = 0; i < mesh.number_of_volume_elements(); i++) 
	{
		Mesh::Volume_element& ve = mesh.find_volume_element(i);
		ofile << std::setw(12) 
			<< ve.ve_conectivity[0] + 1 << " "
			<< ve.ve_conectivity[1] + 1 << " "
			<< ve.ve_conectivity[2] + 1 << " " 
			<< std::endl;
	}
	ofile.close();
	return 1;
}

int Laplace_fem::initial_bc()
{
	Boundary_condition tbc;
	for (int s_id = 0; s_id<mesh.number_of_surface_elements(); ++ s_id)
	{
		Mesh::Surface_element& se = mesh.find_surface_element(s_id);
		if (se.se_tag == 1)
		{
			for (int i = 0; i<se.se_conectivity.size();++i)
			{
				tbc.type = 1;
				tbc.index = se.se_conectivity[0];
				tbc.value = 0;
				bc.push_back(tbc);
			}		
		}
		else if(se.se_tag == 2)
		{
			tbc.type = 2;
			tbc.index = s_id;
			tbc.value = 0;
			bc.push_back(tbc);
		}
	}
	return 0;
}


