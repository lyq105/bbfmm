#include "mesh.h"
#include "quadtree.h"

#include <string>
/// This is an example for build a quad-tree.

void build_a_quadtree( const std::string& fname )
{
	std::string tecfname = fname + ".plt";
	std::string treefname = fname + "_qtree.plt";
	Mesh mesh2d;
	mesh2d.load_mesh( fname.c_str() );
	mesh2d.plot_mesh( tecfname.c_str() );

	typedef qtree_namespace::Quadtree Quadtree;
	Quadtree qtree;
	qtree_namespace::quadtree_creat(qtree,mesh2d);
	qtree_namespace::plot_quadtree(qtree,treefname.c_str());
	qtree_namespace::draw_tree_graph(qtree);
}

int main(int argc, const char *argv[])
{
	std::string fname = argv[1];
	build_a_quadtree(fname);
	return 0;
}
