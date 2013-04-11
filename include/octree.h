//=============================================================================
//   Project name:  head of octree.cpp
//       Filename:  octree.h
//   Descriptions:  四叉树的头文件
//         Author:  liyiqiang(lyq105@163.com)
//         Create:  2011-06-22 21:16:26
//      Last edit:  2011-06-22 21:16:26
//=============================================================================

#ifndef __QUADTREE_H__ 
#define __QUADTREE_H__ 

#include "mesh.h"


namespace octree_namespace{
	typedef Mesh::Point Point;
	// data type's Defination
	typedef struct _treenode
	{
		int isLeaf;   // if the node is a leaf;
		int level;      
		int coord;      // numbering of this node;
		int father;
		int number;
		int startIndex; // start index of element number in the element list
		int maxElem;    // max number of boundary elements in this tree node; 
		double length;
		Point center;
	}TreeNode;

	typedef struct _octree
	{
		int treeType;
		int numberTreenode;   // Number of tree node; 
		int numElem;          // Number of element;
		int treeDepth;        // Depth of tree;
		int* levelCell;    // position of first level l cell; 
		int* numlevelCell; // how many non empty cell of level l;
		int* elemList;    // boundary elements in each tree level, max tree level is 20; 
		TreeNode** treeNodeList; 
	}Octree;


	int octree_creat(Octree& otree, Mesh& mesh);

	int octree_creat_childs(Mesh mesh,Octree& otree,TreeNode* ftnode);
	int octree_init(Octree& otree, int numelem);
	int octree_insert(Octree& otree, TreeNode* tnode);
	int octree_destory(Octree& otree);
	int cal_convex_hall(Mesh& mesh,Point& center,double& length);
	int plot_octree(Octree otree,char filename[]);
	int print_octree_info(Octree otree,char filename[]);

}; // octree_namespace
#endif
