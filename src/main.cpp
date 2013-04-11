#include "lapace2d.h"
#include "timer.h"
#include "quadtree.h"
#include <iostream>
using namespace std;

int test_bem()
{
	Timer tm;
	tm.start();
	Laplace2d laplace("F:\\LiYiQiang\\Desktop\\Fast Bem\\code\\cbem_new\\bem2d.dat");
	tm.stop();
	tm.diplayelapsedtime();
	laplace.initial_bc();
	tm.start();
	laplace.solve();
//	laplace.solve_lapack();
	tm.stop();	
	tm.diplayelapsedtime();
	laplace.print_solution("g:\\bem2d.dat.res.plt");
	return 1;
}

int test_load_mesh() 
{
	Mesh mesh2d,mesh3d;
	mesh2d.load_mesh(std::string("F:\\LiYiQiang\\Desktop\\Fast Bem\\code\\cbem_new\\2d.dat"));
	mesh2d.plot_mesh(std::string("F:\\LiYiQiang\\Desktop\\Fast Bem\\code\\cbem_new\\2d.dat") + ".plt");

	//mesh2d.load_mesh(std::string("e:\\bem2d_big.dat"));
	//mesh2d.plot_mesh(std::string("e:\\bem2d_big.dat") + ".plt");

	mesh3d.load_mesh(std::string("F:\\LiYiQiang\\Desktop\\Fast Bem\\code\\cbem_new\\3d.dat"));
	mesh3d.plot_mesh(std::string("F:\\LiYiQiang\\Desktop\\Fast Bem\\code\\cbem_new\\3d.dat") + ".plt");
	return 1;
}

int test_reformat_mesh_2()
{
	Mesh mesh2d;
	mesh2d.load_mesh("g:\\bem2d.dat");
	mesh2d.save("g:\\bem2d.mesh");
	mesh2d.plot_mesh("g:\\bem2d.plt");
	return 0;
}
int test_load_n_mesh_2()
{
	Mesh mesh2d;
	mesh2d.load_neutral_mesh("g:\\cube.mesh");
//	mesh2d.save("g:\\bem2d.mesh");
	mesh2d.plot_mesh("g:\\cube.plt");
	return 0;
}
int test_quadtree()
{
	Mesh mesh2d;
	mesh2d.load_mesh("g:\\bem2d.dat");
	mesh2d.plot_mesh("g:\\bem2d.plt");

	typedef qtree_namespace::Quadtree Quadtree;
	Quadtree qtree;
	qtree_namespace::quadtree_creat(qtree,mesh2d);
	qtree_namespace::plot_quadtree(qtree,"g:\\bem2d_tree.plt");
	qtree_namespace::draw_tree_graph(qtree);
	return 0;
}


int main()
{	
	//test_load_mesh();
	test_bem();
	//test_reformat_mesh_2();
	//test_load_n_mesh_2();
	//test_quadtree();
	
	return 0;
}