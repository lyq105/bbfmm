#ifndef __elements_H_
#define __elements_H_

///@file Description of elements;


#include "mesh.h"
//#include <mkl_lapacke.h>
#include <vector>
typedef Mesh::Point Point;
typedef Mesh::Volume_element Volume_element;
typedef Mesh::Surface_element Surface_element;

struct Gauss_Int{
	int number_of_gauss_points;
	std::vector<Point> p_list;
	std::vector<double> w_list;
};

// 2d surface element.
class Link2nodeSur
{
public:
	Link2nodeSur(){};
	~Link2nodeSur(){};
	int initial(Mesh mesh, Surface_element ve);
	Point to_real(Point ref_pt);
	double shape(int i, Point ref_pt);
	double dshape(int i,Point ref_pt);
	double dshape_dx(int i,Point ref_pt);
	double get_det_jacobian(){return det_map_jacobian;}
	Gauss_Int gauss_int;

private:
	int initial_nodes(Mesh mesh, Surface_element ve);
	int initial_gauss_int();
	double calcualte_jacobian(Point ref_pt);

	std::vector<Point> real_points_list;
	double map_jacobian[2];
	double inv_jacobian[2];
	double det_map_jacobian;
};

//////////////////////////////////////////////////////////////////////////

// 3d surface element.
class Tri3nodeVol
{
public:
	Tri3nodeVol(){};
	~Tri3nodeVol(){};

	double shape(int i, Point ref_pt);
	double dshape(int i,int d_coord, Point ref_pt);
	double dshape_dx(int i,int d_coord,Point rpt);	
	Point to_real(Point ref_pt);
	double get_det_jacobian(){return det_map_jacobian;};
	int initial(Mesh mesh, Volume_element ve);	
	Gauss_Int gauss_int;
private:

	int initial_nodes(Mesh mesh, Volume_element ve);
	int initial_gauss_int();
	double calcualte_jacobian(Point ref_pt);
	
	std::vector<Point> real_points_list;
	double det_map_jacobian;
	double jacobian[2][2];
	double inv_jacobian[2][2];
};


//////////////////////////////////////////////////////////////////////////

/// 3D Volume element descriptions.

class Tet4nodeVol
{
public:
	Tet4nodeVol(){};
	~Tet4nodeVol();

	Point to_real(Point ref_pt);
	double shape(int i, Point ref_pt);
	double dshape(int i,int d_coord, Point ref_pt);
	double dshape_dx(int i,int d_coord,Point rpt);

	double get_det_jacobian(){return det_map_jacobian;}
	Gauss_Int gauss_int;
	int initial(Mesh mesh, Volume_element ve);
private:

	int initial_nodes(Mesh mesh, Volume_element ve);
	int initial_gauss_int();
	double calcualte_jacobian(Point ref_pt);
private:	

	std::vector<Point> real_points_list;

	double det_map_jacobian;
	double jacobian[3][3];
	double inv_jacobian[3][3];
};

// 3d surface element.
class Tri3nodeSur
{
public:
	Tri3nodeSur(){};
	~Tri3nodeSur();

	double shape(int i, Point ref_pt);
	double dshape(int i,int d_coord, Point ref_pt);
	//double dshape_dx(int i,int d_coord, Point ref_pt);
	Point to_real(Point ref_pt);
	double get_det_jacobian(){return det_map_jacobian;}
	Gauss_Int gauss_int;
	int initial(Mesh mesh, Surface_element ve);
private:

	int initial_nodes(Mesh mesh, Surface_element ve);
	int initial_gauss_int();
	double calcualte_jacobian(Point ref_pt);
	std::vector<Point> real_points_list;
	double jacobian[3][2];
	double det_map_jacobian;
};

#endif
