// @file mesh.h

#ifndef __MESH_H
#define __MESH_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
/// The description of a mesh, which in 2D is constructed by triangles or 
/// tetrahedron, and the surface mesh is constructed by 2d segments or 
/// 3d triangles.

class Mesh
{
public:
	Mesh(){};
	Mesh(std::string filename){load_mesh(filename);};
	~Mesh(){};
	struct Point {
		Point(){ x[0] = x[1] = x[2] = .0;p_tag = 0;}
		Point(const Point& pt){p_tag = pt.p_tag;x[0]=pt.x[0];x[1]=pt.x[1];x[2]=pt.x[2];}
		Point operator=(const Point& pt)
		{ 
			Point rpt;
			rpt.p_tag = pt.p_tag;
			rpt.x[0]=pt.x[0];
			rpt.x[1]=pt.x[1];
			rpt.x[2]=pt.x[2];
			return rpt;
		}
		int p_tag;
		double x[3];
	};
	struct Surface_element {
		int se_mat;
		std::vector<int> se_tag;
		std::vector<int> se_conectivity;
		std::vector<double> se_outer_norm;
	};
	struct Volume_element {
		int ve_mat;
		std::vector<int> ve_tag;
		std::vector<int> ve_conectivity;
	};
	struct Element {
		int e_mat;
		std::vector<int> e_tags;
		std::vector<int> e_conectivity;
	};
	bool load_mesh(std::string filename);
	bool load_gmsh_ascii(std::string filename);
	bool load_neutral_mesh(std::string filename);
	bool save(std::string filename);
	bool print_mesh_info( std::string filename);
	bool plot_mesh(std::string filename);
	int  number_of_point(){return len_point_list;} 
	int  number_of_surface_elements(){return len_sur_elem_list;}
	int  number_of_volume_elements(){return len_vol_elem_list;}

	Point find_point(int p_index){return point_list[p_index];}
	Surface_element find_surface_element(int se_index){return sur_elem_list[se_index];}
	Volume_element find_volume_element(int ve_index){return vol_elem_list[ve_index];}
	Surface_element find_surface_element_by_point(){}

private:
	int mesh_dim;
	int mesh_type;
	int number_of_s_node;
	int number_of_v_node;
	int len_point_list;
	int len_sur_elem_list;
	int len_vol_elem_list;
	std::vector < Point > point_list;
	std::vector < Surface_element > sur_elem_list;
	std::vector < Volume_element > vol_elem_list;
	std::vector < Element > elem_list;

	int cal_surface_outer_norm();	
	int cal_surface_outer_norm_gmsh();
	// Change vertical vector point out of domain.
	int change_to_out_norm(Surface_element& se);
	// Judge is a surface belong to a volume.
	int is_in_volume_element(Volume_element ve, Surface_element se);
	// Find if a volume Point in a surface element.
	int find_surface_element_point(Volume_element ve, Surface_element se);
	int find_parent_volume_element(Surface_element se); 

}; // end of class mesh


#endif // __MESH_H
