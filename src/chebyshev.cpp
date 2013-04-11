

// chebyshev basis function

inline double chebyshev_basis_1d(double x, int num,int Nnum)
{
	double temp=0;
	double omiga=0;
	double x_xi=x+cos((num)/(Nnum-1.0)*PI);//x-xi

	for (int j=0;j<Nnum;j++)
	{
		//omiga_j		
		omiga=pow(-1.0,j);
		if (j==0||j==Nnum-1)
		{
			omiga*=0.5;
		}
		double xj=-cos((j)/(Nnum-1.0)*PI);
		if (j!=num){
			temp+=omiga*x_xi/(x-xj);
		}
		else{
			temp+=omiga;
		}
	}
	//omiga_i
	omiga=pow(-1.0,num);
	if (num==0||num==Nnum-1)
	{
		omiga*=0.5;
	}

	return omiga/temp;
}
// num: 0,...,Nnum*Nnum - 1;
double chebyshev_basis_2d(double x,double y, int num,int Nnum)
{
	int i,j;  // i,j start form 0,not 1. num = j*Nnum + i;
	i=num%Nnum;
	j=int(num/Nnum);
	return chebyshev_basis_1d(x,i,Nnum)*chebyshev_basis_1d(y,j,Nnum);
}

// num: 0,...,Nnum*Nnum*Nnum - 1;
double chebyshev_basis_3d(double x,double y, double z, int num,int Nnum)
{
	int i,j,k;  // i,j start form 0,not 1. num = k*Nnum*Nnum + j*Nnum + i;
	i=num%Nnum;
	j=int(num/Nnum) % Nnum;
	k=int(num/(Nnum*Nnum));
	return chebyshev_basis_1d(x,i,Nnum)
		*chebyshev_basis_1d(y,j,Nnum)
		*chebyshev_basis_1d(z,k,Nnum);
}


// Maps 

Point map_to_ref_1d(Point real_pt, Point center, double length)
{
	Point ref_pt;
	ref_pt.x[0] = 2.0 / length * (real_pt.x[0] - center.x[0]);
	return ref_pt;
}
Point map_to_ref_2d(Point real_pt, Point center, double length)
{
	Point ref_pt;
	ref_pt.x[0] = 2.0 / length * (real_pt.x[0] - center.x[0]);
	ref_pt.x[1] = 2.0 / length * (real_pt.x[1] - center.x[1]);
	return ref_pt;
}
Point map_to_ref_3d(Point real_pt, Point center, double length)
{
	Point ref_pt;
	ref_pt.x[0] = 2.0 / length * (real_pt.x[0] - center.x[0]);
	ref_pt.x[1] = 2.0 / length * (real_pt.x[1] - center.x[1]);
	ref_pt.x[2] = 2.0 / length * (real_pt.x[2] - center.x[2]);
	return ref_pt;
}
Point map_to_real_1d(Point ref_pt, Point center, double length)
{
	Point real_pt;
	real_pt.x[0] = 0.5* length * ref_pt.x[0] + center.x[0];
	return real_pt;
}
Point map_to_real_2d(Point ref_pt, Point center, double length)
{
	Point real_pt;
	real_pt.x[0] = 0.5* length * ref_pt.x[0] + center.x[0];
	real_pt.x[1] = 0.5* length * ref_pt.x[1] + center.x[1];
	return real_pt;
}
Point map_to_real_3d(Point ref_pt, Point center, double length)
{
	Point real_pt;
	real_pt.x[0] = 0.5* length * ref_pt.x[0] + center.x[0];
	real_pt.x[1] = 0.5* length * ref_pt.x[1] + center.x[1];
	real_pt.x[2] = 0.5* length * ref_pt.x[2] + center.x[2];
	return real_pt;
}

typedef double (*Func)(Point pt);

double interpolate_func_1d(Point pt, Func func, Point center, double length, int interpolation_order)
{

	const Point& r_pt = map_to_ref_1d(pt,center,length);

	for (int i = 0; i < interpolation_order; i++)	
	{
		Point ref_pt;
		ref_pt.x[0]=-cos((i)/(interpolation_order-1.0)*PI);

		const Point& real_pt = map_to_real_1d(ref_pt,center,length);
		ret += func(real_pt)*chebyshev_basis_1d(r_pt.x[0] ,r_pt.x[1],i, interpolation_order); 
	}	
	return ret;	
}


double interpolate_func_2d(Point pt, Func func, Point center, double length, int interpolation_order)
{

	const Point& r_pt = map_to_ref_2d(pt,center,length);

	for (int i = 0; i < interpolation_order; i++)	
	{
		for (int j = 0; j < interpolation_order; j++)
		{
			Point ref_pt;
			ref_pt.x[0]=-cos((i)/(interpolation_order-1.0)*PI);
			ref_pt.x[1]=-cos((j)/(interpolation_order-1.0)*PI);

			int num = i*interpolation_order + j;
			const Point& real_pt = map_to_real_2d(ref_pt,center,length);
			ret += func(real_pt)*chebyshev_basis_2d(r_pt.x[0] ,r_pt.x[1],num, interpolation_order); 
		}
	}	
	return ret;	
}


double interpolate_func_3d(Point pt, Func func, Point center, double length, int interpolation_order)
{

	const Point& r_pt = map_to_ref_3d(pt,center,length);

	for (int i = 0; i < interpolation_order; i++)	
	{
		for (int j = 0; j < interpolation_order; j++)
		{
			for (int k = 0; k < interpolation_order; k++)
			{
				Point ref_pt;
				ref_pt.x[0]=-cos((i)/(interpolation_order-1.0)*PI);
				ref_pt.x[1]=-cos((j)/(interpolation_order-1.0)*PI);
				ref_pt.x[2]=-cos((k)/(interpolation_order-1.0)*PI);

				int num = k*interpolation_order*interpolation_order + j*interpolation_order + i;
				const Point& real_pt = map_to_real_3d(ref_pt,center,length);
				ret += func(real_pt)*chebyshev_basis_3d(r_pt.x[0] ,r_pt.x[1],num,	interpolation_order); 
			}
		}	
	}
	return ret;	
}



double inverse_interpolate_func_2d(Point pt, Func func, Point center, double length, int interpolation_order)
{
	return 0;
}

// 上行遍历树结构 计算每一个非空树节点的多极矩
// NOTE:   这里的参数data就是迭代的向量，或者说是矩阵乘向量的向量
//

int quadtree_upward(Quadtree& qtree, Mesh mesh, double* data)
{
	int maxLevel = qtree.treeDepth - 1;  // retrive max tree level.

	// 逐层向根方向遍历树结构 从树的最后一层到第二层
	for (int depth = maxLevel; depth >= 2; --depth)  
	{
		// 获取非空树节点的下标

		int startIndex = qtree.levelCell[ depth ];
		int endIndex = startIndex + qtree.numlevelCell[ depth ] - 1;

		// 遍历同层次的树节点
		for (int cellIndex = startIndex; cellIndex <= endIndex; ++cellIndex) 
		{
			// 这里需要查询父节点的编号  
			// 在传递过程中查询，并传递到父节点
			//{
				//int fatherIndex = qtree.treeNodeList[cellIndex] -> father;
			//}

			// 如果树节点不是叶子节点，则直接使用M2M将其累加到父节点中；
			// 如果该树节点是叶子节点，则计算多极矩，并使用M2M将其累加到其父节点中；
			if ( qtree.treeNodeList[cellIndex] -> isLeaf == 1)  // 是叶子节点
			{
				// 计算多极矩，（调用一个计算接口）， 
				// 输出是该节点的多极矩，输入是该节点的单元编号，
				// 计算编号为cellindex的树节点的多极矩
				// 计算多极矩的时候需要将迭代向量代入计算

				cal_multipole_moment(qtree,cellIndex,mesh,data); 
			}

			//if ( qtree.treeNodeList[cellIndex].isLeaf == 0)
			// 将多极矩传递到父节点中去， {不论当前节点是否是叶子节点}   
			{
				// 传递多极系数使用（M2M公式将当前节点的
				// 多极矩进行转换并累加到父节点中），输入是当前节点的多极矩，

				transfer_mm_to_its_father_by_m2m(qtree,cellIndex);
			}
		} // 遍历层结束

	} // 遍历树结束

	return 0;
}

// 下行遍历树结构，计算局部矩，这里的局部矩包括交互单元和远程单元
// 对整个积分的贡献。

int quadtree_downward(Quadtree& qtree)
{
	// 计算第二层的树节点的多极矩, 使用M2L将交互节点的多极矩转换为局部矩，
	// 累加到中心节点中, 初始计算参数
	{
		// 这一步可以统一到整个树的计算过程中。
	}	
	// 逐层向根方向遍历树结构
	int maxLevel = qtree.treeDepth - 1;  // retrive max tree level.
	for (int depth = 2; depth <= maxLevel; ++ depth)
	{
		int startIndex = qtree.levelCell[ depth ];
		int endIndex = startIndex + qtree.numlevelCell[ depth ] - 1;

		// 遍历同层次的树节点
		for (int cellIndex = startIndex; cellIndex <= endIndex; ++cellIndex) 
		{
			QuadtreeNode& qnode = *qtree.treeNodeList[cellIndex];
			// 遍历交互节点列表，使用M2L，将局部矩累加到节点上
			// 需要查询交互节点
			for (int i = 0; i < qnode.lenInterList; i++)
			{
				int iIndex = qnode.interList[i];
				transfer_lm_by_m2l(qtree, cellIndex, iIndex);	
			}
			// 将远程单元的贡献累加到当前单元中，
			// 只要使用L2L将父节点的局部系数转移到当前节点的中心上
			{
				if (depth == 2 ) continue; // 如果当前节点是第二层节点，则不必使用L2L
				transfer_lm_by_l2l(qtree, cellIndex);
			}
		} // 遍历层

	} // 遍历树
}

