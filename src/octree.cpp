//=============================================================================
//   Project name:  libtree
//       Filename:  octree.cpp
//   Descriptions:  constrction a quad tree code.
//         Author:  liyiqiang(lyq105@163.com)
//         Create:  2011-06-22 21:16:26
//      Last edit:  2011-06-22 21:16:26
//=============================================================================

#include <stdio.h>
#include <stdlib.h>

#include "mesh.h"
#include "octree.h"
#include <iostream>
using namespace std;


namespace octree_namespace{

	static int maxTreeDepth = 50;
	static int maxEleminCell = 2;
	static int main_number = 0;

	// 计算网格点的一个立方体凸包
	int cal_convex_hall(Mesh& mesh,Point& center,double& length)
	{
		int i,j;
		double x_max[3],x_min[3];
		length = 0;
		for (j = 0; j < 3; j++)
		{
			x_max[j] = x_min[j]=0;
			for (i = 0; i < mesh.number_of_point(); i++)
			{
				Mesh::Point& pt = mesh.find_point(i);
				x_max[j] = x_max[j] >= pt.x[j]? x_max[j] : pt.x[j];
				x_min[j] = x_min[j] <= pt.x[j]? x_min[j] : pt.x[j];
			}
			center.x[j] = 0.5*(x_max[j] + x_min[j]);
			length = length >= (x_max[j] - x_min[j]) ? length: (x_max[j] - x_min[j]);
		}
		length *= 1.1;
		return 0;
	}

	int octree_creat(Octree& otree, Mesh& mesh)
	{
		int i;
		int depth, nodenumber;
		int startIndex, endIndex;

		TreeNode* root = (TreeNode*)malloc(sizeof(TreeNode));

		cal_convex_hall(mesh,root->center,root->length);


		root -> number = 0;    // root 's global numbering
		root -> level = 0;     // root 's level
		root -> isLeaf = 0;    // is a leaf node
		root -> father = -1;   // root father's numbering is -1 
		root -> startIndex = 0; 
		root -> maxElem = mesh.number_of_surface_elements();
		root -> coord = 0;

		// 初始化树
		octree_init(otree,root->maxElem);
		// 插入根结点
		octree_insert(otree,root);

		otree.levelCell[0] = 0;
		otree.numlevelCell[0] = 1;

		if (root -> maxElem <= maxEleminCell)
		{
			return 1;
		}

		// 逐层遍历生成树结构

		for (depth = 1; depth < maxTreeDepth; depth++)
		{
			nodenumber = otree.numberTreenode;
			// 遍历depth-1层的树节点
			startIndex = otree.levelCell[depth - 1] ; // depth-1层节点开始的下标
			endIndex = startIndex + otree.numlevelCell[depth - 1] - 1; // depth-1层节点结束的下标

			for (i = startIndex; i <= endIndex; i++)
			{
				if (otree.treeNodeList[i] -> maxElem > maxEleminCell)
				{
					otree.treeNodeList[i] -> isLeaf = 0;
					octree_creat_childs(mesh,otree,otree.treeNodeList[i]); // 四分节点，并将非空树节点插入树节点列表
				}
				else
					otree.treeNodeList[i] -> isLeaf = 1;
			}

			if (otree.numberTreenode != nodenumber)	 // 如果有新节点产生,
			{
				otree.levelCell[depth] = endIndex + 1;
				otree.numlevelCell[depth] = otree.numberTreenode - nodenumber;
			}
			else break;
		}
		otree.treeDepth = depth;

		print_octree_info(otree,"treeinfo.txt");
		plot_octree(otree,"treeinfo.plt");

		return 0;
	}

	int octree_creat_childs(Mesh mesh,Octree& otree,TreeNode* ftnode)
	{
		int i,j;
		int pos[8];
		int startIndex,endIndex,eindex,estartIndex;

		// 分配一个和父节点单元数一样的临时矩阵
		int** temparray = new int*[8];
		for (i = 0; i < 8; i++) {
			temparray[i] = new int[ftnode -> maxElem];
		}

		double coef_x[8]={-0.25,0.25,0.25,-0.25,-0.25,0.25,0.25,-0.25};
		double coef_y[8]={-0.25,-0.25,0.25,0.25,-0.25,-0.25,0.25,0.25};
		double coef_z[8]={-0.25,-0.25,-0.25,-0.25,0.25,0.25,0.25,0.25};

		// 计算每一个象限的单元数
		for (i = 0; i < 8; i++) {
			pos[i] = 0;
		}

		// startIndex和endIndex表示在树的单元列表中开始和结束索引。
		startIndex = ftnode -> startIndex;
		endIndex = startIndex + ftnode -> maxElem;


		for(i = startIndex; i< endIndex; i++)  //  遍历父节点中的单元
		{
			// 计算单元中心的坐标
			eindex = otree.elemList[i];

			double x[3];
			for (int cindex = 0; cindex < 3; cindex++) 
			{
				for (int pindex = 0; pindex < 3; pindex++) 
				{
					Mesh::Surface_element& se = mesh.find_surface_element(eindex);
					Mesh::Point& pt1 = mesh.find_point(se.se_conectivity[0]);
					x[cindex] += pt1.x[cindex];
				}
				x[cindex] *= 0.5;
			}

			// 计算8个区域中的单元数
			if(x[2] < ftnode -> center.x[2])
			{ // U
				if( x[1] >= ftnode -> center.x[1])
				{ // N
					if (x[0] >= ftnode -> center.x[0])
					{ // E
						temparray[2][pos[2]] = eindex; // 将该单元编号存入临时矩阵
						pos[2]++;
					}
					else
					{  // W
						temparray[3][pos[3]] = eindex;
						pos[3]++;
					}
				}
				else
				{	// S
					if (x[0] >= ftnode -> center.x[0])
					{ // E
						temparray[1][pos[1]] = eindex; // 将该单元编号存入临时矩阵
						pos[1]++;
					}
					else
					{ // W
						temparray[0][pos[0]] = eindex; // 将该单元编号存入临时矩阵
						pos[0]++;
					}
				}
			}
			else
			{
				if( x[2] >= ftnode -> center.x[1])
				{ // N
					if (x[0] >= ftnode -> center.x[0])
					{ // E
						temparray[6][pos[6]] = eindex; // 将该单元编号存入临时矩阵
						pos[6]++;
					}
					else
					{  // W
						temparray[7][pos[7]] = eindex;
						pos[7]++;
					}
				}
				else
				{	// S
					if (x[0] >= ftnode -> center.x[0])
					{ // E
						temparray[5][pos[5]] = eindex; // 将该单元编号存入临时矩阵
						pos[5]++;
					}
					else
					{ // W
						temparray[4][pos[4]] = eindex; // 将该单元编号存入临时矩阵
						pos[4]++;
					}
				}
			}

		}

		estartIndex = ftnode->startIndex;
		for (j = 0; j < 8; j++)
		{
			if (pos[j] != 0 ) // 创建非空节点
			{
				TreeNode* childnode = new TreeNode;
				// 赋层数；
				childnode -> level = ftnode -> level + 1;
				childnode -> coord = ftnode -> coord * 8 + j;
				// 父节点
				childnode -> father = ftnode -> number;
				// 节点的整体编号；
				main_number += 1;
				childnode -> number = main_number;
				// 单元开始下标；
				childnode -> startIndex = estartIndex;
				estartIndex += pos[j];
				// 节点中的单元个数；
				childnode -> maxElem = pos[j];
				// 单元中心点；
				childnode -> center.x[0] = ftnode -> center.x[0] + coef_x[j]* ftnode -> length;
				childnode -> center.x[1] = ftnode -> center.x[1] + coef_y[j]* ftnode -> length;
				childnode -> center.x[2] = ftnode -> center.x[1] + coef_z[j]* ftnode -> length;

				// 单元边长；
				childnode -> length = 0.5 * ftnode -> length;

				for (int k = childnode->startIndex;k<childnode->startIndex + pos[j]; k++)
				{
					otree.elemList[k] = temparray[j][k - childnode->startIndex];
				}
				// 将节点插入节点列表；
				octree_insert(otree,childnode);
			}
		}

		for (i = 0; i < 8; i++)
		{
			delete [] temparray[i];
		}
		delete [] temparray;
		return 0;
	}
	int octree_init(Octree& otree, int numelem)
	{
		// 初始化一个空树
		otree.numberTreenode = 0;
		otree.numElem = numelem;
		otree.treeDepth = 0;

		otree.levelCell = new int[maxTreeDepth];
		otree.numlevelCell = new int[maxTreeDepth];

		otree.elemList = new int[numelem];

		for (int i = 0; i < numelem; i++) {
			otree.elemList[i] = i;
		}

		otree.treeNodeList = NULL;
		return 0;
	}

	// 将节点插入节点列表
	int octree_insert(Octree& otree, TreeNode* tnode)
	{
		int i;
		otree.numberTreenode += 1;
		TreeNode** templist;
		templist = new TreeNode*[otree.numberTreenode];

		for (i = 0; i < otree.numberTreenode - 1; i++) {
			templist[i] = otree.treeNodeList[i];
		}
		templist[otree.numberTreenode - 1] = tnode;

		// 删除旧的树节点列表
		delete [] otree.treeNodeList;

		otree.treeNodeList = templist;
		return 0;
	}

	// 释放树占据的空间
	int octree_destory(Octree& otree)
	{
		int i;
		delete [] otree.levelCell;
		delete [] otree.numlevelCell;
		delete [] otree.elemList;

		for (i = 0; i < otree.numberTreenode; i++) {
			//delete [] otree.treeNodeList[i]->elementList;
			delete otree.treeNodeList[i];
		}
		delete [] otree.treeNodeList;

		return 0;
	}

	int print_octree_info(Octree otree,char filename[])
	{
		FILE* fp;
		fp = fopen(filename,"w");
		// Export tree nodes;
		fprintf(fp,"number isleaf  coord    father   level    center_x\tcenter_y\tlength\tmaxnum\n");
		for(int i=0; i<otree.numberTreenode; i++)
		{
			fprintf(fp,"%8d %8d %8d %8d %8d %12.8f %12.8f %12.8f %d\n",otree.treeNodeList[i]->number,otree.treeNodeList[i]->isLeaf,otree.treeNodeList[i]->coord,otree.treeNodeList[i]->father,otree.treeNodeList[i]->level,otree.treeNodeList[i]->center.x[0],otree.treeNodeList[i]->center.x[1],
				otree.treeNodeList[i]->length,otree.treeNodeList[i]->maxElem);
		}
		fclose(fp);
		return 1;
	}
	int plot_octree(Octree otree,char filename[])
	{
		FILE* fp;
		double coef_x[8]={-0.25,0.25,0.25,-0.25,-0.25,0.25,0.25,-0.25};
		double coef_y[8]={-0.25,-0.25,0.25,0.25,-0.25,-0.25,0.25,0.25};
		double coef_z[8]={-0.25,-0.25,-0.25,-0.25,0.25,0.25,0.25,0.25};

		fp =fopen(filename,"w");
		fprintf(fp,"zone N=%d,E=%d,T=\"octree\",F=FEPOINT,ET=cube\n",4*otree.numberTreenode,otree.numberTreenode);
		for(int i=0; i<otree.numberTreenode; i++)
		{
			for(int j =0; j<8;j++)
			{
				double x = otree.treeNodeList[i]->center.x[0] + coef_x[j] * otree.treeNodeList[i]->length;
				double y = otree.treeNodeList[i]->center.x[1] + coef_y[j] * otree.treeNodeList[i]->length;
				double z = otree.treeNodeList[i]->center.x[1] + coef_z[j] * otree.treeNodeList[i]->length;
				fprintf(fp,"%12.8f %12.8f\n",x,y,z);
			}
		}
		for(int i=0; i<otree.numberTreenode; i++)
		{
			fprintf(fp, "%d %d %d %d %d %d %d %d\n",8*i+1,8*i+2,8*i+3,8*i+4,8*i+5,8*i+6,8*i+7,8*i+8);
		}
		return 0;
	}

}; // namespace octree