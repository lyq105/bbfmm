#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <math.h>
#include "mesh.h"


int Mesh::cal_surface_outer_norm()
{
	if (mesh_dim == 2)  // for 2d triangle element.
	{
		double nx = 0;
		double ny = 0;
		double norm_length;

		///  A vertical vector is like 
		//  \f[
		//  	\frac{1}{\sqrt{a_1^2+a_2^2}}(-a_2,a_1).
		//  \f]
		for (int i = 0; i < len_sur_elem_list; i++)
		{
			Surface_element& se = sur_elem_list[i];		
			Point& n0 = point_list[se.se_conectivity[0]]; 
			Point& n1 = point_list[se.se_conectivity[1]]; 

			nx = - (n0.x[1] - n1.x[1]);
			ny =   (n0.x[0] - n1.x[0]);
			norm_length = sqrt(nx*nx + ny*ny);
			se.se_outer_norm.push_back( nx / norm_length);
			se.se_outer_norm.push_back( ny / norm_length);	
			change_to_out_norm(se);
		}
	}

	if (mesh_dim == 3)
	{
		double nx = 0;
		double ny = 0;
		double nz = 0;
		double norm_length;

		for (int i = 0; i < len_sur_elem_list; i++)
		{
			Surface_element& se = sur_elem_list[i];		
			Point& n0 = point_list[se.se_conectivity[0]]; 
			Point& n1 = point_list[se.se_conectivity[1]]; 
			Point& n2 = point_list[se.se_conectivity[2]];

			Point p1,p2,p3;
			p1.x[0] = n1.x[0] - n0.x[0];
			p1.x[1] = n1.x[1] - n0.x[1];
			p1.x[2] = n1.x[2] - n0.x[2];
			p2.x[0] = n2.x[0] - n0.x[0];
			p2.x[1] = n2.x[1] - n0.x[1];
			p2.x[2] = n2.x[2] - n0.x[2];

			// cross product of p1 and p2.

			p3.x[0] =    p1.x[1] * p2.x[2] - p1.x[2] * p2.x[1];
			p3.x[1] = - (p1.x[0] * p2.x[2] - p1.x[2] * p2.x[0]);
			p3.x[2] =    p1.x[0] * p2.x[1] - p1.x[1] * p2.x[0];

			norm_length = sqrt(p3.x[0]*p3.x[0] + p3.x[1]*p3.x[1] + p3.x[2]*p3.x[2]);
			se.se_outer_norm.push_back( p3.x[0] / norm_length);
			se.se_outer_norm.push_back( p3.x[1] / norm_length);
			se.se_outer_norm.push_back( p3.x[2] / norm_length);	

			change_to_out_norm(se);
		}
	}
	return 1;
}
int Mesh::cal_surface_outer_norm_gmsh()
{
	if (mesh_dim == 2)  // for 2d triangle element.
	{
		double nx = 0;
		double ny = 0;
		double norm_length;

		///  A vertical vector is like 
		//  \f[
		//  	\frac{1}{\sqrt{a_1^2+a_2^2}}(-a_2,a_1).
		//  \f]
		for (int i = 0; i < len_sur_elem_list; i++)
		{
			Surface_element& se = sur_elem_list[i];		
			Point& n0 = point_list[se.se_conectivity[0]]; 
			Point& n1 = point_list[se.se_conectivity[1]]; 

			nx = - (n0.x[1] - n1.x[1]);
			ny =   (n0.x[0] - n1.x[0]);
			norm_length = sqrt(nx*nx + ny*ny);
			se.se_outer_norm.push_back( nx / norm_length);
			se.se_outer_norm.push_back( ny / norm_length);	
//			change_to_out_norm(se);
		}
	}

	if (mesh_dim == 3) ///TODO:
	{
		double nx = 0;
		double ny = 0;
		double nz = 0;
		double norm_length;

		for (int i = 0; i < len_sur_elem_list; i++)
		{
			Surface_element& se = sur_elem_list[i];		
			Point& n0 = point_list[se.se_conectivity[0]]; 
			Point& n1 = point_list[se.se_conectivity[1]]; 
			Point& n2 = point_list[se.se_conectivity[2]];

			Point p1,p2,p3;
			p1.x[0] = n1.x[0] - n0.x[0];
			p1.x[1] = n1.x[1] - n0.x[1];
			p1.x[2] = n1.x[2] - n0.x[2];
			p2.x[0] = n2.x[0] - n0.x[0];
			p2.x[1] = n2.x[1] - n0.x[1];
			p2.x[2] = n2.x[2] - n0.x[2];

			// cross product of p1 and p2.

			p3.x[0] =    p1.x[1] * p2.x[2] - p1.x[2] * p2.x[1];
			p3.x[1] = - (p1.x[0] * p2.x[2] - p1.x[2] * p2.x[0]);
			p3.x[2] =    p1.x[0] * p2.x[1] - p1.x[1] * p2.x[0];

			norm_length = sqrt(p3.x[0]*p3.x[0] + p3.x[1]*p3.x[1] + p3.x[2]*p3.x[2]);
			se.se_outer_norm.push_back( p3.x[0] / norm_length);
			se.se_outer_norm.push_back( p3.x[1] / norm_length);
			se.se_outer_norm.push_back( p3.x[2] / norm_length);	

//			change_to_out_norm(se);
		}
	}
	return 1;
}
// Change vertical vector point out of domain.
int Mesh::change_to_out_norm(Surface_element& se)
{
	//std::cout << "invokes change norm\n";
	// find out which volume element contains se.
	Volume_element& ve = vol_elem_list[find_parent_volume_element(se)];

	// find out the node which is not in the surface element.
	Point& n3 = point_list[find_surface_element_point(ve,se)];
	//		std::cout<< "ve = " << find_surface_element_Point(ve,se) << std::endl;

	Point& n0 = point_list[se.se_conectivity[0]];
	if (mesh_dim == 2)
	{
		double nx = n3.x[0] - n0.x[0];
		double ny = n3.x[1] - n0.x[1];	
		//			std::cout << n3.x[0] << " === " << n3.x[1] <<std::endl;
		if (nx* se.se_outer_norm[0] + ny * se.se_outer_norm[1] > 0)
		{
			se.se_outer_norm[0] *= -1;
			se.se_outer_norm[1] *= -1;
		}
	}
	if (mesh_dim == 3)
	{
		double nx = n3.x[0] - n0.x[0];
		double ny = n3.x[1] - n0.x[1];
		double nz = n3.x[2] - n0.x[2];	

		if (nx* se.se_outer_norm[0] + ny * se.se_outer_norm[1]+ nz * se.se_outer_norm[2]> 0)
		{
			se.se_outer_norm[0] *= -1;
			se.se_outer_norm[1] *= -1;
			se.se_outer_norm[2] *= -1;
		}
	}
	return 0;
}

// Find surface in which volume element.
int Mesh::find_parent_volume_element(Surface_element se)
{
	for (int i = 0; i< len_vol_elem_list; ++i )
	{
		Volume_element& ve = vol_elem_list[i];
		if (is_in_volume_element(ve,se) == 1)
			return i;  // return the number of volume element.
	}
	return -1;
}

// Judge is a surface belong to a volume.
int Mesh::is_in_volume_element(Volume_element ve, Surface_element se)
{
	int sum = 0;
	for (int i = 0; i < number_of_s_node; i++ )
		for (int j = 0; j < number_of_v_node; j++ )
		{
			//				std::cout << se.se_conectivity[i] << "vs " << ve.ve_conectivity[j] << "\n";
			if (se.se_conectivity[i] == ve.ve_conectivity[j])
			{
				sum ++;	
			}
		}
		if (sum == number_of_s_node) return 1;
		return 0;
}

// Find if a volume Point in a surface element.
int Mesh::find_surface_element_point(Volume_element ve, Surface_element se)
{
	// find node which is not in surface element,but in volume element 
	// which contain this surface element. 
	for (int i = 0; i < number_of_v_node; i++)
	{
		int sum = 0;
		for (int j = 0; j < number_of_s_node; j++)
		{
			// 				std::cout << se.se_conectivity[j] << "vs " << ve.ve_conectivity[i] << "\n";

			if (ve.ve_conectivity[i] != se.se_conectivity[j]){
				sum ++;
			}
		}
		if (sum == number_of_s_node)
		{
			return ve.ve_conectivity[i]; 
		}
	}
	return -1;
}
/// load a mesh from a mesh data file.
bool Mesh::load_mesh(std::string filename)
{
	std::string buffer;
	std::istringstream buffer_parser;
	std::ifstream infile(filename.c_str(),std::ios::in);
	if ( !infile.is_open() )
	{
		std::cout << "ERROR #1: Can not load the mesh data file \""
			<< filename << "\"!!"
			<< std::endl;
		//exit(0);
	}
	std::getline(infile,buffer,'\n');
	if (buffer.length() == 0 || buffer[0] == '#')
	{
		std::getline(infile,buffer,'\n');
	}
	buffer_parser.str(buffer);
	buffer_parser >> mesh_dim 
		>> len_point_list
		>> len_sur_elem_list
		>> len_vol_elem_list;

	if (mesh_dim == 2)  // for 2d triangle element.
	{
		number_of_s_node = 2;
		number_of_v_node = 3;
	}
	else if (mesh_dim == 3)  // for 3d tetrahedron element.
	{
		number_of_s_node = 3;
		number_of_v_node = 4;
	}

	//	std::cout << len_point_list;

	// The next data block is the Point_coordinates_list; 
	for (int i = 0; i < len_point_list; ++i)
	{
		Point ptemp;
		std::getline(infile,buffer,'\n');
		if (buffer.length() == 0 || buffer[0] == '#')
		{
			std::getline(infile,buffer,'\n');
		}
		buffer_parser.clear();
		buffer_parser.str(buffer);
		buffer_parser >> ptemp.p_tag;
		for (int p_id = 0; p_id < mesh_dim; ++ p_id){			
			buffer_parser >> ptemp.x[p_id];
		}
		//		std::cout << ptemp.x[0] << "   " << ptemp.x[1] << std::endl;
		point_list.push_back(ptemp);
	}

	// The next data block is the surface element list; 

	for (int i = 0; i < len_sur_elem_list; i++) 
	{
		Surface_element se;
		std::getline(infile,buffer,'\n');
		if (buffer.length() == 0 || buffer[0] == '#')
		{
			std::getline(infile,buffer,'\n');
		}
		int tag;
		buffer_parser.clear();
		buffer_parser.str(buffer);
		buffer_parser >> se.se_mat >>tag;
		//		std::cout << se.se_tag << "  " << se.se_mat << std::endl;
		for (int se_id = 0; se_id < number_of_s_node; ++ se_id)
		{
			int temp;
			buffer_parser >> temp;
			
			se.se_conectivity.push_back(temp-1);
// 			std::cout << se.se_conectivity[se_id] << " ";
// 			std::cout << temp-1 << " ";
		}
// 		std::cout << std::endl;
		//		std::cout << " length is "<< se.se_conectivity.size() << std::endl;
		sur_elem_list.push_back(se);	
	}
	std::cout << " surface mesh have "<< sur_elem_list.size() << " elements "<< std::endl;
	// The next data block is the volume element list; 
	for (int i = 0; i < len_vol_elem_list; i++) 
	{
		Volume_element ve;
		int tag;
		std::getline(infile,buffer,'\n');
		if (buffer.length() == 0 || buffer[0] == '#')
		{
			std::getline(infile,buffer,'\n');
		}
		buffer_parser.clear();
		buffer_parser.str(buffer);
		buffer_parser >> ve.ve_mat >> tag;
		for (int ve_id = 0; ve_id < number_of_v_node; ++ ve_id)
		{
			int temp;
			buffer_parser >> temp;
			ve.ve_conectivity.push_back(temp-1);
		}
		vol_elem_list.push_back(ve);	
	}
	std::cout << " Volume mesh have "<< vol_elem_list.size() << " elements "<< std::endl;

	cal_surface_outer_norm(); //
	infile.close();
	return 1;
}

bool Mesh::load_gmsh_ascii(std::string filename)
{
	std::string buffer;
	std::istringstream buffer_parser;
	std::ifstream infile(filename.c_str(),std::ios::in);
	
	if ( !infile.is_open() ){
		std::cout << "ERROR #1: Can not load the mesh data file \""
			<< filename << "\"!!"
			<< std::endl;
		return 0;
	}
	// get format info
	while(!infile.eof()){
		std::getline(infile,buffer,'\n');
		std::string start_tag = "$MeshFormat";
		std::string end_tag = "$EndMeshFormat";
		if (buffer == start_tag){
			std::getline(infile,buffer,'\n');
			while (buffer != end_tag){
				std::cout << buffer << std::endl;
				std::getline(infile,buffer,'\n');
			}
			break;
		}
	}
	// get nodes information	
	infile.seekg(0);
	while(!infile.eof()){
		int npoints = 0,nnum;
		Point ptemp;
		std::string start_tag = "$Nodes";
		std::string end_tag = "$EndNodes";
		
		std::getline(infile,buffer,'\n');
		if (buffer == start_tag){
			std::getline(infile,buffer,'\n');
			buffer_parser.clear();
			buffer_parser.str(buffer);
			buffer_parser >> npoints;
			while (buffer != end_tag){
				std::getline(infile,buffer,'\n');
				buffer_parser.clear();
				buffer_parser.str(buffer);
				buffer_parser >> nnum >>ptemp.x[0] >> ptemp.x[1] >> ptemp.x[2];
				point_list.push_back(ptemp);			
			}
			break;
		}
	}
	// get elements information.
	infile.seekg(0);
	while(!infile.eof()){
		int nelems,num,etype,netags,tags;
		std::string start_tag = "$Elements";
		std::string end_tag = "$EndElemsnts";

		std::getline(infile,buffer,'\n');

		if (buffer == start_tag){
			std::getline(infile,buffer,'\n');
			buffer_parser.clear();
			buffer_parser.str(buffer);
			buffer_parser >> nelems;

			while (buffer != end_tag){
				Element elem;
				std::getline(infile,buffer,'\n');
				buffer_parser.clear();
				buffer_parser.str(buffer);
				buffer_parser >> num >> etype >> netags;

				for (int i=0;i< netags;i++){
					int etag;
					buffer_parser >> etag;
					elem.e_tags.push_back(etag);
				}
				
				while(!buffer_parser.eof())
				{
					int econ;
					buffer_parser >> econ;
					elem.e_conectivity.push_back(econ);
				}
				elem_list.push_back(elem);
			}
			break;
		}
	}
	
	infile.close();
	return 1;
}


// Load mesh from Netgen neutral mesh format.
bool Mesh::load_neutral_mesh(std::string filename)
{
	std::string buffer;
	std::istringstream buffer_parser;
	std::ifstream infile(filename.c_str(),std::ios::in);
	if ( !infile.is_open() )
	{
		std::cout << "ERROR #1: Can not load the mesh data file \""
			<< filename << "\"!!"
			<< std::endl;
		//exit(0);
	}

	mesh_dim = 3;
	number_of_s_node = 3;
	number_of_v_node = 4;

	std::getline(infile,buffer,'\n');
	if (buffer.length() == 0 || buffer[0] == '#')
	{
		std::getline(infile,buffer,'\n');
	}
	buffer_parser.str(buffer);
	buffer_parser >> len_point_list;


	// The next data block is the Point_coordinates_list; 
	for (int i = 0; i < len_point_list; ++i)
	{
		Point ptemp;
		std::getline(infile,buffer,'\n');
		if (buffer.length() == 0 || buffer[0] == '#')
		{
			std::getline(infile,buffer,'\n');
		}
		buffer_parser.clear();
		buffer_parser.str(buffer);

		for (int p_id = 0; p_id < mesh_dim; ++ p_id){			
			buffer_parser >> ptemp.x[p_id];
		}
		//		std::cout << ptemp.x[0] << "   " << ptemp.x[1] << std::endl;
		point_list.push_back(ptemp);
	}

	// The next data block is the volume element list; 
	std::getline(infile,buffer,'\n');
	if (buffer.length() == 0 || buffer[0] == '#')
	{
		std::getline(infile,buffer,'\n');
	}
	buffer_parser.clear();
	buffer_parser.str(buffer);
	buffer_parser >> len_vol_elem_list;
	//std::cout << len_vol_elem_list;
	for (int i = 0; i < len_vol_elem_list; i++) 
	{
		Volume_element ve;

		std::getline(infile,buffer,'\n');
		if (buffer.length() == 0 || buffer[0] == '#')
		{
			std::getline(infile,buffer,'\n');
		}
		buffer_parser.clear();
		buffer_parser.str(buffer);
		int ve_tag=0;
		buffer_parser >> ve_tag ;
		ve.ve_tag.push_back(ve_tag);
		for (int ve_id = 0; ve_id < number_of_v_node; ++ ve_id)
		{
			int temp;
			buffer_parser >> temp;
			ve.ve_conectivity.push_back(temp-1);
		}
		vol_elem_list.push_back(ve);	
	}
	std::cout << " Volume mesh have "<< vol_elem_list.size() << " elements "<< std::endl;

	// The next data block is the surface element list; 
	std::getline(infile,buffer,'\n');
	if (buffer.length() == 0 || buffer[0] == '#')
	{
		std::getline(infile,buffer,'\n');
	}
	buffer_parser.clear();
	buffer_parser.str(buffer);
	buffer_parser >> len_sur_elem_list;

	for (int i = 0; i < len_sur_elem_list; i++) 
	{
		Surface_element se;
		std::getline(infile,buffer,'\n');
		if (buffer.length() == 0 || buffer[0] == '#')
		{
			std::getline(infile,buffer,'\n');
		}

		buffer_parser.clear();
		buffer_parser.str(buffer);
		int se_tag=0;
		buffer_parser >> se_tag;
		se.se_tag.push_back(se_tag);
		//		std::cout << se.se_tag << "  " << se.se_mat << std::endl;
		for (int se_id = 0; se_id < number_of_s_node; ++ se_id)
		{
			int temp;
			buffer_parser >> temp;
			se.se_conectivity.push_back(temp-1);
			// 			std::cout << se.se_conectivity[se_id] << " ";
			// 			std::cout << temp-1 << " ";
		}
		// 		std::cout << std::endl;
		//		std::cout << " length is "<< se.se_conectivity.size() << std::endl;
		sur_elem_list.push_back(se);	
	}
	std::cout << " surface mesh have "<< sur_elem_list.size() << " elements "<< std::endl;


	cal_surface_outer_norm_gmsh(); //
	infile.close();
	return 1;
}
bool Mesh::print_mesh_info( std::string filename )
{
	std::ofstream ofile(filename.c_str(),std::ios::out);
	std::cout << "This mesh contains " << len_point_list << "points\n";
	ofile << "This mesh contains " << len_point_list << "points\n";
	return 1;
}

bool Mesh::plot_mesh(std::string filename)
{
	std::ofstream ofile(filename.c_str(),std::ios::out);
// Plot volume mesh.

	ofile << "zone f=fepoint,e=" 
		<< number_of_volume_elements()
		<< ", n= " 
		<< number_of_point();
	if (mesh_dim == 3)
	{
		ofile<< ",et=tetrahedron"<< std::endl;
	}
	if (mesh_dim == 2)
	{
		ofile<< ",et=triangle"<< std::endl;
	}		

	for (int i = 0; i< number_of_point(); i++)
	{
		for (int p_id = 0; p_id < mesh_dim; ++ p_id){
			ofile << std::setw(12) << std::setprecision(8)
				<< find_point(i).x[p_id];
		}
		ofile << std::endl;
	}
	for (int i = 0; i < number_of_volume_elements(); i++) 
	{
		for (int ve_id = 0; ve_id < number_of_v_node; ++ ve_id)
		{
			ofile << std::setw(12) 
				<< find_volume_element(i).ve_conectivity[ve_id] + 1; 
		}
		ofile << std::endl;
	}

// Plot surface mesh.

	ofile << "zone f=fepoint,e=" 
		<< number_of_surface_elements()
		<< ", n= " 
		<< number_of_point()
		<< ",et=triangle"
		<< std::endl;


	for (int i = 0; i< number_of_point(); i++)
	{
		for (int p_id = 0; p_id < mesh_dim; ++ p_id){
			ofile << std::setw(12) << std::setprecision(8)
				  << find_point(i).x[p_id];
		}
		ofile << std::endl;
	}
	for (int i = 0; i < number_of_surface_elements(); i++) 
	{
		for (int se_id = 0; se_id < number_of_s_node; ++ se_id){
			ofile << std::setw(12) 
				<< find_surface_element(i).se_conectivity[se_id] + 1; 
//			std::cout << se.se_conectivity[se_id] + 1 ;
		}
// 		std::cout << std::endl;
		if (mesh_dim == 2){
			ofile << std::setw(12) 
				    << find_surface_element(i).se_conectivity[0] + 1;
		}
		ofile << std::endl;
	}
	ofile.close();

	std::cout << "Mesh is write to " << filename
		<< "\n Please use tecplot to open it. "
		<< std::endl;
	return 0;
}

// Save mesh data in Neutral Format
bool Mesh::save( std::string filename )
{
	std::ofstream ofile(filename.c_str(),std::ios::out);
	ofile << this -> number_of_point() << std::endl;
	for (int i = 0; i < this->number_of_point();i++ )
	{
		for (int j =0; j< this->mesh_dim; ++j)
		{
			ofile << this->point_list[i].x[j] << "   ";
		}
		ofile << std::endl;
	}
	ofile << this -> number_of_volume_elements() << std::endl;
	for (int i = 0; i < this->number_of_volume_elements(); i++) 
	{
		for (int ve_id = 0; ve_id < number_of_v_node; ++ ve_id)
		{
			ofile << std::setw(12) 
				<< find_volume_element(i).ve_conectivity[ve_id] + 1; 
		}
		ofile << std::endl;
	}
	ofile << this -> number_of_surface_elements() << std::endl;
	for (int i = 0; i < this->number_of_surface_elements();i++ )
	{
		for (int se_id = 0; se_id < number_of_s_node; ++ se_id){
			ofile << std::setw(12) 
				<< find_surface_element(i).se_conectivity[se_id] + 1; 
			//			std::cout << se.se_conectivity[se_id] + 1 ;
		}
		// 		std::cout << std::endl;
		if (mesh_dim == 2)
			//ofile << std::setw(12) << se.se_conectivity[0] + 1;
			ofile << std::endl;
	}
	return 1;
}
/*
int Mesh::load_mesh_ansys(std::string& nodefilename,std::string elemfilename)
{
	FILE* fp;
	FILE* efp;

	char buffer[200];
	char x[3][21]={'\0'};
	char elem[8][10]={'\0'};
	char ch;

	fp=fopen(nodefilename,"r");
	if (fp == NULL)
	{
		printf("Cannot open file %s!!\n",nodefilename);
		return 0;
	}
	efp=fopen(elemfilename,"r");
	if (efp == NULL)
	{
		printf("Cannot open file %s!!\n",elemfilename);
		return 0;
	}

	int nnum=0;
	while((ch=fgetc(fp))!=EOF)
	{
		if(ch == '\n')nnum++;
	}
	rewind(fp);


	int elemnum=0;
	while((ch=fgetc(efp))!=EOF)
	{
		if(ch=='\n')elemnum++;
	}
	rewind(efp);

	number_of_s_node = nnum;
	mesh->TotleElementNumber = elemnum;

#ifdef _USE_CPP_
	mesh->MeshPointList = new Point[mesh->TotleNodeNumber];
	mesh->CellList = new Cell[mesh->TotleElementNumber];
#else  // using Ansi C.
	mesh->MeshPointList = (Point*)malloc(sizeof(Point)*mesh->TotleNodeNumber);
	mesh->CellList = (Point*)malloc(sizeof(Cell)*mesh->TotleElementNumber);
#endif


	int pindex = 0;
	while(fgets(buffer,200,fp)!=NULL)
	{
		for (int j=0;j<mesh->dim;j++)
		{
			x[j][0]='\0';
			for(int i=0;i<20;i++)
			{
				int strid=j*20+i+8;
				if(buffer[strid] != '\0')
				{
					x[j][i] = buffer[strid];
					x[j][i+1]='\0';
				}
			}
			mesh->MeshPointList[pindex].x[j]=atof(x[j]);
		}

		pindex++;

//#ifdef _DEBUG_
		//system("pause");
		printf("%8d ",pindex-1);
		for (int j=0;j<mesh->dim;j++)
			printf("%12.8f ",atof(x[j]));
		printf("\n");
//#endif
	}

	fclose(fp);

	int eindex=0;
	while(fgets(buffer,200,efp)!=NULL)
	{

		// Allocate memory of every cell.

#ifdef _USE_CPP_
		mesh->CellList[eindex].cellnumber = new int[mesh->nodepercell];
#else  // using Ansi C.
		mesh->CellList[eindex].cellnumber = (int*) malloc(sizeof(int)*mesh->nodepercell);
#endif
		for (int j=0;j < mesh->nodepercell;j++)
		{
			elem[j][0]='\0';
			for(int i=0;i<efmt;i++)
			{
				int strid=j*efmt+i;
				if(et == tet3d && j == 3)
				{
					strid = 4*efmt+i;
				}
				if(buffer[strid] != '\0')
				{
					elem[j][i] = buffer[strid];
					elem[j][i+1]='\0';
				}
			}
			mesh->CellList[eindex].cellnumber[j]=atoi(elem[j]);
		}

		eindex++; // element index ++;

#ifdef _DEBUG_
		for (int j=0;j< nodenum;j++)
			printf("%s ",(elem[j]));
		printf("\n");
#endif
	}

	fclose(efp);
	return 1;
}
*/
// int main()
// {	
// 	Mesh mesh2d,mesh3d;
// // 	mesh2d.load_mesh(std::string("F:\\LiYiQiang\\Desktop\\Fast Bem\\code\\cbem_new\\2d.dat"));
// // 	mesh2d.plot_mesh(std::string("F:\\LiYiQiang\\Desktop\\Fast Bem\\code\\cbem_new\\2d.dat") + ".plt");
// 	mesh2d.load_mesh(std::string("e:\\bem2d_big.dat"));
// 	mesh2d.plot_mesh(std::string("e:\\bem2d_big.dat") + ".plt");
// // 	mesh3d.load_mesh(std::string("F:\\LiYiQiang\\Desktop\\Fast Bem\\code\\cbem_new\\3d.dat"));
// // 	mesh3d.plot_mesh(std::string("F:\\LiYiQiang\\Desktop\\Fast Bem\\code\\cbem_new\\3d.dat") + ".plt");
// 
// 	return 0;
// }
