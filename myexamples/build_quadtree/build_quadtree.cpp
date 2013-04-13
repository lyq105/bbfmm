#include "mesh.h"
#include "quadtree.h"


/// This is an example for build a quad-tree.

void build_a_quadtree()
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

int main(int argc, const char *argv[])
{
	build_a_quadtree();
	return 0;
}
