#include "elements.h"


#include "mesh.h"
//#include <mkl_lapacke.h>
#include <vector>

template <typename _T>
double inverse_matrix(const _T& omat, _T& imat,int n)
{
	double** mat = new double*[n];

	for(int i = 0; i<n ; ++i)
	{
		mat[i] = new double[n];
		for(int j=0;j<n;j++)
		{
			mat[i][j]= omat[i][j];
			imat[i][j] = 0.0;
		}
		imat[i][i] = 1.0;
	}

	for(int i = 0; i<n-1; ++i)
	{
		for (int k=1; k<n-i;k++)
		{
			double factor = mat[i+k][i]/mat[i][i];
			for (int j = 0;j<n;++j)
			{
				mat[i+k][j] -= factor*mat[i][j]; 
				imat[i+k][j] -= factor*imat[i][j]; 
			}
		}
	}

	for(int i = n-1; i>=0; --i)
	{
		for (int k=1; k<=i;++k)
		{
			double factor = mat[i-k][i]/mat[i][i];
			for (int j = 0;j<n;++j)
			{
				mat[i-k][j] -= factor*mat[i][j]; 
				imat[i-k][j] -= factor*imat[i][j]; 
			}
		}
	}
	double det = 1;
	for (int i=0;i<n;i++){
		det*= mat[i][i];
		for(int j=0;j<n;j++){
			imat[i][j] /= mat[i][i];
		}
	}
	return det;
}

double Link2nodeSur::shape(int i, Point ref_pt)
{
	if(i == 0) return 0.5 * (1 + ref_pt.x[0]);
	if(i == 1) return 0.5 * (1 - ref_pt.x[0]);
	return 0;
}

double Link2nodeSur::dshape(int i,Point ref_pt)
{
	if(i == 0) return  0.5;
	if(i == 1) return -0.5;
	return 0;
}
double Link2nodeSur::dshape_dx(int i,Point rpt)
{
	return dshape(i,rpt)*inv_jacobian[i]; 
}

Point Link2nodeSur::to_real(Point ref_pt)
{
	Point real_pt;
	for (int s_id = 0; s_id < 2; ++s_id)
	{
		real_pt.x[0] += real_points_list[s_id].x[0] * shape(s_id,ref_pt);
		real_pt.x[1] += real_points_list[s_id].x[1] * shape(s_id,ref_pt);
	}
	return real_pt;
};

int Link2nodeSur::initial_nodes( Mesh mesh, Surface_element se )
{
	int point_index;
	for (int p_id = 0; p_id < 2; ++p_id)
	{
		point_index = se.se_conectivity[p_id];
		real_points_list.push_back(mesh.find_point(point_index));
	}
	return 0;
}

int Link2nodeSur::initial_gauss_int()
{
	gauss_int.number_of_gauss_points = 2;
	Point p1,p2;
	p1.x[0] = -1.0/sqrt(3.0);
	p2.x[0] =  1.0/sqrt(3.0); 
	gauss_int.p_list.push_back(p1);
	gauss_int.p_list.push_back(p2);
	gauss_int.w_list.push_back(1.0);
	gauss_int.w_list.push_back(1.0);

	return 0;
}


//////////////////////////////////////////////////////////////////////////
//    \frac{dx} {d\yita} \frac{ dy} { d\yita}
//////////////////////////////////////////////////////////////////////////

double Link2nodeSur::calcualte_jacobian(Point ref_pt)
{
	double h = 1e-6;
	double det_jocabian;
	Point ref_pt1,ref_pt2;
	ref_pt1.x[0] = ref_pt.x[0] + h;
	ref_pt2.x[1] = ref_pt.x[1] + h;

	map_jacobian[0] = (to_real(ref_pt1).x[0] - to_real(ref_pt).x[0])/h;
	map_jacobian[1] = (to_real(ref_pt1).x[1] - to_real(ref_pt).x[1])/h;
	inv_jacobian[0] = 1.0/map_jacobian[0];
	inv_jacobian[1] = 1.0/map_jacobian[1];

	det_jocabian = sqrt(map_jacobian[0] * map_jacobian[0] + map_jacobian[1] * map_jacobian[1]);
	return det_jocabian;
}

int Link2nodeSur::initial( Mesh mesh, Surface_element ve )
{
	initial_gauss_int();
	initial_nodes(mesh,ve);
	Point ref_pt;
	calcualte_jacobian(ref_pt);
	return 0;
}

//////////////////////////////////////////////////////////////////////////

double Tri3nodeVol::shape(int i, Point ref_pt)
{
	if(i == 0) return ref_pt.x[0];
	if(i == 1) return ref_pt.x[1];
	if(i == 2) return 1 - ref_pt.x[0] - ref_pt.x[1];
	return 0;
}

double Tri3nodeVol::dshape(int i,int d_coord, Point ref_pt)
{
	if(i == 0){
		if (d_coord == 0) return 1;
		if (d_coord == 1) return 0;
	}
	if(i == 1){
		if (d_coord == 0) return 0;
		if (d_coord == 1) return 1;
	}
	if(i == 2) {
		if (d_coord == 0) return -1;
		if (d_coord == 1) return -1;
	}
	return 0;
}
double Tri3nodeVol::dshape_dx(int i,int d_coord,Point rpt)
{
	return dshape(i,0,rpt)*inv_jacobian[0][d_coord] 
	+ dshape(i,1,rpt)*inv_jacobian[1][d_coord];
}
int Tri3nodeVol::initial_nodes(Mesh mesh, Volume_element ve)
{
	int point_index;
	for (int p_id = 0; p_id < 3; ++p_id)
	{
		point_index = ve.ve_conectivity[p_id];
		real_points_list.push_back(mesh.find_point(point_index));
	}
	return 0;
};

/*
\f[
\frac{\partial x_i}{\partial \eta_j} 
\f]
because of the shape function is a polynomial function, so we can use numerical difference 
to calculate jacobian.
*/


Point Tri3nodeVol::to_real(Point ref_pt)
{
	Point real_pt;
	for (int s_id = 0; s_id < 3; ++s_id)
	{
		real_pt.x[0] += real_points_list[s_id].x[0] * shape(s_id,ref_pt);
		real_pt.x[1] += real_points_list[s_id].x[1] * shape(s_id,ref_pt);
	}
	return real_pt;
};

double Tri3nodeVol::calcualte_jacobian( Point ref_pt )
{
	double h = 1e-6;
	double det_jocabian;
	Point ref_pt1,ref_pt2;
	ref_pt1.x[0] = ref_pt.x[0] + h;
	ref_pt1.x[1] = ref_pt.x[1];
	ref_pt2.x[0] = ref_pt.x[0];
	ref_pt2.x[1] = ref_pt.x[1] + h;

	jacobian[0][0] = (to_real(ref_pt1).x[0] - to_real(ref_pt).x[0])/h;
	jacobian[0][1] = (to_real(ref_pt1).x[1] - to_real(ref_pt).x[1])/h;
	jacobian[1][0] = (to_real(ref_pt1).x[0] - to_real(ref_pt).x[0])/h;
	jacobian[1][1] = (to_real(ref_pt1).x[1] - to_real(ref_pt).x[1])/h;

	det_jocabian = jacobian[0][0] *	jacobian[1][1] - jacobian[0][1] *	jacobian[1][0];

	inv_jacobian[0][0] = jacobian[1][1]/det_jocabian;
	inv_jacobian[0][1] = - jacobian[0][1]/det_jocabian;
	inv_jacobian[1][1] = jacobian[0][0]/det_jocabian;
	inv_jacobian[1][0] = - jacobian[1][0]/det_jocabian;

	return det_jocabian;
}

int Tri3nodeVol::initial_gauss_int()
{
	double gauss_pt_x[3] ={0.0,0.5,0.5}; 
	double gauss_pt_y[3] ={0.5,0.5,0.0}; 
	double gauss_wt[3] ={1/3.0,1/3.0,1/3.0};

	gauss_int.number_of_gauss_points = 3;
	for (int i=0;i<3;++i)
	{
		Point pt;
		pt.x[0] = gauss_pt_x[i];
		pt.x[1] = gauss_pt_y[i];
		gauss_int.p_list.push_back(pt);
		gauss_int.w_list.push_back(gauss_wt[i]);
	}
	return 0;
}
int Tri3nodeVol::initial( Mesh mesh, Volume_element ve )
{
	initial_gauss_int();
	initial_nodes(mesh,ve);
	Point ref_pt;
	calcualte_jacobian(ref_pt);
	return 0;
}

//////////////////////////////////////////////////////////////////////////


double Tet4nodeVol::shape(int i, Point ref_pt)
{
	if(i == 0) return ref_pt.x[0];
	if(i == 1) return ref_pt.x[1];
	if(i == 2) return ref_pt.x[2];
	if(i == 3) return 1 - ref_pt.x[0] - ref_pt.x[1] - ref_pt.x[2];
	return 0;
}

double Tet4nodeVol::dshape(int i,int d_coord, Point ref_pt)
{
	if(i == 0){
		if (d_coord == 0) return 1;
		if (d_coord == 1) return 0;
		if (d_coord == 2) return 0;
	}
	if(i == 1){
		if (d_coord == 0) return 0;
		if (d_coord == 1) return 1;
		if (d_coord == 2) return 0;
	}
	if(i == 2) {
		if (d_coord == 0) return 0;
		if (d_coord == 1) return 0;
		if (d_coord == 2) return 1;
	}
	if(i == 3) {
		if (d_coord == 0) return -1;
		if (d_coord == 1) return -1;
		if (d_coord == 2) return -1;
	}
	return 0;
}

double Tet4nodeVol::dshape_dx(int i,int d_coord,Point rpt)
{
	return dshape(i,0,rpt)*inv_jacobian[0][d_coord] 
	     + dshape(i,1,rpt)*inv_jacobian[1][d_coord]
	     + dshape(i,2,rpt)*inv_jacobian[2][d_coord];
}

int Tet4nodeVol::initial_nodes(Mesh mesh, Volume_element ve)
{
	int point_index;
	for (int p_id = 0; p_id < 4; ++ p_id)
	{
		point_index = ve.ve_conectivity[p_id];
		real_points_list.push_back(mesh.find_point(point_index));
	}
	return 0;
};
/*
\f[
\frac{\partial x_i}{\partial \eta_j} 
\f]
Because of the shape function is a polynomial function, so we can use numerical difference 
to calculate jacobian.
*/

double Tet4nodeVol::calcualte_jacobian(Point ref_pt)
{

	for (int i = 0; i < 3; ++ i)
	{
		for (int j = 0; j < 3; ++ j)
		{
			jacobian[i][j] = 0;
			for (int k=0; k < 4; ++k)
			{
				jacobian[i][j] += real_points_list[k].x[j]*dshape(k,j,ref_pt);
			}
		}
	}
	//det_jocabian = determinate(3,(double**)jacobian);
// 求逆矩阵
  det_map_jacobian = inverse_matrix<double[3][3]>(jacobian,inv_jacobian,3);
	
	return det_map_jacobian;
};

Point Tet4nodeVol::to_real(Point ref_pt)
{
	Point real_pt;
	for (int s_id = 0; s_id < 3; ++s_id)
	{
		real_pt.x[0] += real_points_list[s_id].x[0] * shape(s_id,ref_pt);
		real_pt.x[1] += real_points_list[s_id].x[1] * shape(s_id,ref_pt);
		real_pt.x[2] += real_points_list[s_id].x[2] * shape(s_id,ref_pt);
	}
	return real_pt;
};

int Tet4nodeVol::initial_gauss_int()
{
	double a = 0.50541020;
	double b = 0.13319680;
	double gauss_pt_x[4] ={a,b,b,b}; 
	double gauss_pt_y[4] ={b,a,b,b}; 
	double gauss_pt_z[4] ={b,b,a,b};
	double gauss_wt[4] ={1/4.0, 1/4.0, 1/4.0,1/4.0};

	gauss_int.number_of_gauss_points = 4;
	for (int i=0;i<4;++i)
	{
		Point pt;
		pt.x[0] = gauss_pt_x[i];
		pt.x[1] = gauss_pt_y[i];
		pt.x[2] = gauss_pt_z[i];
		gauss_int.p_list.push_back(pt);
		gauss_int.w_list.push_back(gauss_wt[i]);
	}
	return 0;
}
int Tet4nodeVol::initial( Mesh mesh, Volume_element ve )
{
	initial_gauss_int();
	initial_nodes(mesh,ve);
	Point ref_pt;
	calcualte_jacobian(ref_pt);
	return 0;
}

//////////////////////////////////////////////////////////////////////////
double Tri3nodeSur::shape(int i, Point ref_pt)
{
	if(i == 0) return ref_pt.x[0];
	if(i == 1) return ref_pt.x[1];
	if(i == 2) return 1 - ref_pt.x[0] - ref_pt.x[1];
	return 0;
}

double Tri3nodeSur::dshape(int i,int d_coord, Point ref_pt)
{
	if(i == 0){
		if (d_coord == 0) return 1;
		if (d_coord == 1) return 0;
	}
	if(i == 1){
		if (d_coord == 0) return 0;
		if (d_coord == 1) return 1;
	}
	if(i == 2) {
		if (d_coord == 0) return -1;
		if (d_coord == 1) return -1;
	}
	return 0;
}
//double Tri3nodeSur::dshape_dx(int i,int d_coord, Point ref_pt)
//{
	//double ret=0;
	//for (int k = 0; k<3; k++)
	//{
		 //ret += dshape(i,k,ref_pt)*inv_jacobian[k][d_coord];
	//}
	//return ret;
//}

int Tri3nodeSur::initial_nodes(Mesh mesh, Surface_element ve)
{
	int point_index;
	for (int p_id = 0; p_id < 3; ++p_id)
	{
		point_index = ve.se_conectivity[p_id];
		real_points_list.push_back(mesh.find_point(point_index));
	}
	return 0;
};
/*
\f[
\frac{\partial x_i}{\partial \eta_j} 
\f]
Because of the shape function is a polynomial function, so we can use numerical difference 
to calculate jacobian.
*/

double Tri3nodeSur::calcualte_jacobian(Point ref_pt)
{
	for (int i = 0; i < 3; ++i){
		for (int j = 0; j < 2; ++j){
			for (int k = 0; k < 3; ++k){
				jacobian[i][j] += real_points_list[k].x[j]*dshape(k,j,ref_pt);
			}
		}
	}
	double d1,d2,d3;
	double d1_mat[2][2];
	double d2_mat[2][2];
	double d3_mat[2][2];

	d1_mat[0][0] = jacobian[1][0];d1_mat[0][1] = jacobian[1][1];
	d1_mat[1][0] = jacobian[2][0];d1_mat[1][1] = jacobian[2][1];

	d2_mat[0][0] = jacobian[0][0];d2_mat[0][1] = jacobian[0][1];
	d2_mat[1][0] = jacobian[1][0];d2_mat[1][1] = jacobian[1][1];

	d3_mat[0][0] = jacobian[0][0];d3_mat[0][1] = jacobian[0][1];
	d3_mat[1][0] = jacobian[2][0];d3_mat[1][1] = jacobian[2][1];

	d1 = d1_mat[0][0]*d1_mat[1][1] - d1_mat[1][0]*d1_mat[0][1];
	d2 = d2_mat[0][0]*d2_mat[1][1] - d2_mat[1][0]*d2_mat[0][1];
	d3 = d3_mat[0][0]*d3_mat[1][1] - d3_mat[1][0]*d3_mat[0][1];

	return sqrt(d1*d1 + d2*d2 + d3*d3);
};

Point Tri3nodeSur::to_real(Point ref_pt)
{
	Point real_pt;
	for (int s_id = 0; s_id < 3; ++s_id)
	{
		real_pt.x[0] += real_points_list[s_id].x[0] * shape(s_id,ref_pt);
		real_pt.x[1] += real_points_list[s_id].x[1] * shape(s_id,ref_pt);
		real_pt.x[2] += real_points_list[s_id].x[1] * shape(s_id,ref_pt);
	}
	return real_pt;
};

int Tri3nodeSur::initial_gauss_int()
{
	double gauss_pt_x[3] ={0.0,0.5,0.5}; 
	double gauss_pt_y[3] ={0.5,0.5,0.0}; 
	double gauss_wt[3] ={1/3.0,1/3.0,1/3.0};
	
	gauss_int.number_of_gauss_points = 3;
	for (int i=0;i<3;++i)
	{
		Point pt;
		pt.x[0] = gauss_pt_x[i];
		pt.x[1] = gauss_pt_y[i];
		gauss_int.p_list.push_back(pt);
		gauss_int.w_list.push_back(gauss_wt[i]);
	}
	return 0;
}
int Tri3nodeSur::initial( Mesh mesh, Surface_element se )
{
	initial_gauss_int();
	initial_nodes(mesh,se);
	Point ref_pt;
	calcualte_jacobian(ref_pt);
	return 0;
}
