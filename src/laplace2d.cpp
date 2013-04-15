#include "lapace2d.h"
#include "matrix.h"
#include <iomanip>
//#include <mkl_lapacke.h>
// define some variables

void GAUSSELIMINATE(double **E,int M,double *RHS);

// describe the boundary condition
int Laplace2d::initial_bc()
{
	bc.freedom = bmesh.number_of_surface_elements();
	//std::cout << bc.freedom << std::endl;
 	bc.bctags.resize(bc.freedom);
 	bc.bcvalues.resize(bc.freedom);
//int test = 0 ;
	for (int i = 0; i < bc.freedom; i++)
	{	
		const Mesh::Surface_element& se = bmesh.find_surface_element(i);
		if(se.se_tag[0] == 1)    // first kind boundary condition
		{
// 			test++;
			bc.bctags[i] = 1;
			bc.bcvalues[i] = 200;
		}
		if(se.se_tag[0] == 2) // second kind boundary condition
		{ 
// 			test++;
			bc.bctags[i] = 2;
			bc.bcvalues[i] = 10;
		}
	}

//	for (int i = 0; i< bmesh.number_of_point(); ++i)
//	{
//		std::cout 
//			<< bmesh.find_point(i).p_tag << " " 
//			<< bmesh.find_point(i).x[0]  << " "
//			<< bmesh.find_point(i).x[1]  << " "
//			<< std::endl;
//	}
//	std::cout << "test = " << test << std::endl;
	return 1;
}

int Laplace2d::solve()
{
	int i,j;
	double hij,gij,x,y;

	double ** matrix;
	double * rhs;

	const double pi = 3.14159265;
	double gaussnum = 4;
	double gausspt[4] = { 0.86113631, -0.86113631, 0.33998014, -0.33998014 };
	double gausswt[4] = {0.34785485,0.34785485,0.65214515,0.65214515};

	//double gaussnum = 7;
	//double gausspt[7] = {0, 0.94910791, -0.94910791, 0.74153119, -0.74153119, 0.40584515,-0.40584515 };
	//double gausswt[7] = {0.41795918,0.12948497,0.12948497,0.27970539,0.27970539,0.38183005,0.38183005};

	int freedom = bmesh.number_of_surface_elements();

	matrix = new double* [freedom];
	for (int i = 0; i < freedom; i++) 
	{
		matrix[i]= new double[freedom];
		for (int j = 0; j < freedom; j++) {
			matrix[i][j] = 0;
		}
	}

	rhs = new double[freedom];

	for (i = 0; i < freedom; i++) 
	{
		rhs[i] = 0;
	}

	for(i = 0; i< freedom; i++) 
	{
		const Mesh::Surface_element& source_e = bmesh.find_surface_element(i);
		const Mesh::Point& s1 = bmesh.find_point(source_e.se_conectivity[0]);
		const Mesh::Point& s2 = bmesh.find_point(source_e.se_conectivity[1]);
		// find out the source point
		double source_x0 = 0.5 * (s1.x[0] + s2.x[0]);
		double source_x1 = 0.5 * (s1.x[1] + s2.x[1]);
//		printf("%f %f\n",source_x0,source_x1);
		for (j = 0; j < freedom; j++) 
		{
			// find out the number of 2 ends of a cell.
			const Mesh::Surface_element& se = bmesh.find_surface_element(j);
			const Mesh::Point& p1 = bmesh.find_point(se.se_conectivity[0]);
			const Mesh::Point& p2 = bmesh.find_point(se.se_conectivity[1]);
//			std::cout << se.se_conectivity[0]

			// calculate the length of a cell.
			double gama_j = sqrt((p1.x[0] - p2.x[0])*(p1.x[0] - p2.x[0]) 
				+ (p1.x[1] - p2.x[1])*(p1.x[1] - p2.x[1]));
			// calculate the outer unit normal vector of the cell.
			double unx =  se.se_outer_norm[0];
			double uny =  se.se_outer_norm[1];

//			printf("%f %f\n",unx,uny);

			if(i == j) // Diagonal element of matrix which have singular integration.
			{
				hij = 1.0/2;
				gij = gama_j / (2*pi) * (log(2.0/gama_j) + 1);
			}
			else
			{
				// use gauss integral to calculate
				// \f[
				//  {h_{ij}} = \int_{{\Gamma _e}}^{} {F(x,y)} d{\Gamma _e} 
				//     = \sum\limits_{}^{} {F(x,{y_i})*0.5*d*{w_i}}
				//  {g_{ij}} = \int_{{\Gamma _e}}^{} {G(x,y)} d{\Gamma _e}
				//    = \sum\limits_{}^{} {G(x,{y_i})*0.5*d*{w_i}}
				// \f]
				hij = 0;
				gij = 0;
				for(int gindex = 0; gindex< gaussnum; gindex++)
				{
					x = 0.5*(p1.x[0] - p2.x[0])*gausspt[gindex] + 0.5*(p1.x[0] + p2.x[0]);
					y = 0.5*(p1.x[1] - p2.x[1])*gausspt[gindex] + 0.5*(p1.x[1] + p2.x[1]);
					double rr = ((x-source_x0)*(x-source_x0) + (y-source_x1)*(y-source_x1));
					hij += 1.0/ rr * gausswt[gindex];
					gij += - log(sqrt(rr)) *gausswt[gindex];
				}
				double d = (x-source_x0)*unx + (y-source_x1)*uny;
				hij *= -1.0/(4*pi)*gama_j*d;
				gij *= 1.0/(4*pi)*gama_j;
			}
			// Assemble the matrix and right hand side vector;
			if(bc.bctags[j] == 1)
			{
				matrix[i][j] += -gij;
				rhs[i] += -hij*bc.bcvalues[j];
			}
			else if (bc.bctags[j] == 2)
			{
				matrix[i][j] += hij;
				rhs[i] += gij*bc.bcvalues[j];
			}
		}
	} // loop for source point .
	
	GAUSSELIMINATE(matrix,freedom,rhs);

	solution.q.resize(freedom,0);
	solution.u.resize(freedom,0);
	
	for (i = 0; i < freedom; i++) 
	{
//		std::cout << rhs[i]<<std::endl;
		if (bc.bctags[i] == 1)
		{
			solution.u[i] = bc.bcvalues[i];
			solution.q[i] = rhs[i];
		}
		if (bc.bctags[i] == 2)
		{
			solution.u[i] = rhs[i];
//			std::cout << "u == " <<solution.u[i] << std::endl;
			solution.q[i] = bc.bcvalues[i];
		}
//		std::cout << solution.u[i] << " " << solution.q[i] << std::endl;
	}

	//TODO: free memory of matrix and rhs.
	delete [] rhs;
	rhs = NULL;
	for (int i = 0; i< freedom; ++i)
	{
		delete [] matrix[i];
	}
	delete [] matrix;
	matrix = NULL;
	return 0;
};

int Laplace2d::solve_lapack()
{
	int i,j;
	double hij,gij,x,y;

	double * matrix;
	double * rhs;

	const double pi = 3.14159265;
	double gaussnum = 4;
	double gausspt[4] = { 0.86113631, -0.86113631, 0.33998014, -0.33998014 };
	double gausswt[4] = {0.34785485,0.34785485,0.65214515,0.65214515};

	//double gaussnum = 7;
	//double gausspt[7] = {0, 0.94910791, -0.94910791, 0.74153119, -0.74153119, 0.40584515,-0.40584515 };
	//double gausswt[7] = {0.41795918,0.12948497,0.12948497,0.27970539,0.27970539,0.38183005,0.38183005};

	int freedom = bmesh.number_of_surface_elements();

	matrix = new double [freedom*freedom];
	

	rhs = new double[freedom];

	for (i = 0; i < freedom; i++) 
	{
		rhs[i] = 0;
	}

	for(i = 0; i< freedom; i++) 
	{
		const Mesh::Surface_element& source_e = bmesh.find_surface_element(i);
		const Mesh::Point& s1 = bmesh.find_point(source_e.se_conectivity[0]);
		const Mesh::Point& s2 = bmesh.find_point(source_e.se_conectivity[1]);
		// find out the source point
		double source_x0 = 0.5 * (s1.x[0] + s2.x[0]);
		double source_x1 = 0.5 * (s1.x[1] + s2.x[1]);

		for (j = 0; j < freedom; j++) 
		{
			// find out the number of 2 ends of a cell.
			const Mesh::Surface_element& se = bmesh.find_surface_element(j);
			const Mesh::Point& p1 = bmesh.find_point(se.se_conectivity[0]);
			const Mesh::Point& p2 = bmesh.find_point(se.se_conectivity[1]);
			//			std::cout << se.se_conectivity[0]

			// calculate the length of a cell.
			double gama_j = sqrt((p1.x[0] - p2.x[0])*(p1.x[0] - p2.x[0]) 
				+ (p1.x[1] - p2.x[1])*(p1.x[1] - p2.x[1]));
			// calculate the outer unit normal vector of the cell.
			double unx =  se.se_outer_norm[0];
			double uny =  se.se_outer_norm[1];

			if(i == j) // Diagonal element of matrix which have singular integration.
			{
				hij = 1.0/2;
				gij = gama_j / (2*pi) * (log(2.0/gama_j) + 1);
			}
			else
			{
				// use gauss integral to calculate
				// \f[
				//  {h_{ij}} = \int_{{\Gamma _e}}^{} {F(x,y)} d{\Gamma _e} 
				//     = \sum\limits_{}^{} {F(x,{y_i})*0.5*d*{w_i}}
				//  {g_{ij}} = \int_{{\Gamma _e}}^{} {G(x,y)} d{\Gamma _e}
				//    = \sum\limits_{}^{} {G(x,{y_i})*0.5*d*{w_i}}
				// \f]
				hij = 0;
				gij = 0;
				for(int gindex = 0; gindex< gaussnum; gindex++)
				{
					x = 0.5*(p1.x[0] - p2.x[0])*gausspt[gindex] + 0.5*(p1.x[0] + p2.x[0]);
					y = 0.5*(p1.x[1] - p2.x[1])*gausspt[gindex] + 0.5*(p1.x[1] + p2.x[1]);
					double rr = ((x-source_x0)*(x-source_x0) + (y-source_x1)*(y-source_x1));
					hij += 1.0/ rr * gausswt[gindex];
					gij += - log((rr)) *gausswt[gindex];
				}
				double d = (x-source_x0)*unx + (y-source_x1)*uny;
				hij *= -1.0/(4*pi)*gama_j*d;
				gij *= 1.0/(4*pi)*gama_j;
			}
			// Assemble the matrix and right hand side vector;
			if(bc.bctags[j] == 1)
			{
				matrix[i + j *freedom] += -gij;
				rhs[i] += -hij*bc.bcvalues[j];
			}
			else if (bc.bctags[j] == 2)
			{
				matrix[i + j *freedom] += hij;
				rhs[i] += gij*bc.bcvalues[j];
			}
		}
	} // loop for source point .
	// use Intel MKL to solve this linear equation.

	int* ipiv = new int[freedom];
	int n = freedom;
	int nrhs = 1;
	int lda = freedom;
	int ldb = freedom;

//	int info = LAPACKE_dgesv(LAPACK_COL_MAJOR,n,nrhs,matrix,lda,ipiv,rhs,ldb);
	delete [] ipiv;

//	if (info == 0){
//		std::cout << " Successfully call Lapack_dgesv!!"  << std::endl;
//	}
//	else{
//		std::cout << "Calling lapack_dgesv failed!!"  << std::endl;
//	}
//	solution.q.resize(freedom,0);
//	solution.u.resize(freedom,0);

	for (i = 0; i < freedom; i++) 
	{
		//		std::cout << rhs[i]<<std::endl;
		if (bc.bctags[i] == 1)
		{
			solution.u[i] = bc.bcvalues[i];
			solution.q[i] = rhs[i];
		}
		if (bc.bctags[i] == 2)
		{
			solution.u[i] = rhs[i];
			//			std::cout << "u == " <<solution.u[i] << std::endl;
			solution.q[i] = bc.bcvalues[i];
		}
		//		std::cout << solution.u[i] << " " << solution.q[i] << std::endl;
	}

	//TODO: free memory of matrix and rhs.
	delete [] rhs;
	rhs = NULL;
	delete [] matrix;
	matrix = NULL;
	return 0;
};

int Laplace2d::print_solution(std::string resfname)
{
	std::ofstream ofile(resfname.c_str(),std::ios::out);

	ofile << "zone f=fepoint,e=" 
		<< bmesh.number_of_volume_elements()
		<< ", n= " 
		<< bmesh.number_of_point()
		<< ",et=triangle"
		<< std::endl;
		

	for (int i = 0; i< bmesh.number_of_point(); i++)
	{
		const Mesh::Point& p0 = bmesh.find_point(i);
		//std::cout << p0.p_tag << std::endl;
		if(p0.p_tag == 0)   // boundary point
		{
			ofile << std::setw(18) << std::setprecision(8)
				//<< p0.p_tag << " " 
				  << p0.x[0] << " "
			      << p0.x[1] << " "
			      << u_boundary(i)
				  << std::endl;
		}
		else{
			ofile << std::setw(18) << std::setprecision(10)
				  << p0.x[0] << " "
		    	  << p0.x[1] << " "
			      << u_field(p0.x[0],p0.x[1])
				  << std::endl;

		}
	}
	for (int i = 0; i < bmesh.number_of_volume_elements(); i++) 
	{
		const Mesh::Volume_element& ve = bmesh.find_volume_element(i);
		ofile << std::setw(12) 
			  << ve.ve_conectivity[0] + 1 << " "
			  << ve.ve_conectivity[1] + 1 << " "
			  << ve.ve_conectivity[2] + 1 << " " 
			  << std::endl;
	}
	ofile.close();
	return 0;
}

double Laplace2d::u_boundary(int number)
{
	double k = 0;
	double sum = 0;
	int freedom = bmesh.number_of_surface_elements();

	for(int i = 0; i<freedom;i++)
	{
		const Mesh::Surface_element& se = bmesh.find_surface_element(i);
 		if(number == se.se_conectivity[0] || number == se.se_conectivity[1])
		{
			k += 1;
			sum += solution.u[i];
		}
	}
	return sum/k;
}
double Laplace2d::u_field(double sx0, double sx1)
{
	int j;
	double unx,uny;
	double hij,gij;
	double rr,x,y,d;
	double ret = 0;

	const double pi = 3.14159265;
	double gaussnum = 4;
	double gausspt[4] = { 0.86113631, -0.86113631, 0.33998014, -0.33998014 };
	double gausswt[4] = {0.34785485,0.34785485,0.65214515,0.65214515};
// 	double gaussnum = 7;
// 	double gausspt[7] = {0, 0.94910791, -0.94910791, 0.74153119, -0.74153119, 0.40584515,-0.40584515 };
// 	double gausswt[7] = {0.41795918,0.12948497,0.12948497,0.27970539,0.27970539,0.38183005,0.38183005};
	int freedom = bmesh.number_of_surface_elements();

 	for (j = 0; j < freedom; j++) 
 	{
		// find out the number of 2 ends of a cell.
		const Mesh::Surface_element& se = bmesh.find_surface_element(j);
		const Mesh::Point& p1 = bmesh.find_point(se.se_conectivity[0]);
		const Mesh::Point& p2 = bmesh.find_point(se.se_conectivity[1]);

		// calculate the length of a cell.
		double gama_j = sqrt((p1.x[0] - p2.x[0])*(p1.x[0] - p2.x[0]) 
			+ (p1.x[1] - p2.x[1])*(p1.x[1] - p2.x[1]));

		// calculate the outer unit normal vector of the cell.
		unx =  se.se_outer_norm[0];
		uny =  se.se_outer_norm[1];

		hij = 0;
		gij = 0;
		for(int gindex = 0; gindex< gaussnum; gindex++)
		{
			x = 0.5*(p1.x[0]-p2.x[0])*gausspt[gindex] + 0.5*(p1.x[0] + p2.x[0]);
			y = 0.5*(p1.x[1]-p2.x[1])*gausspt[gindex] + 0.5*(p1.x[1] + p2.x[1]);
			rr = ((x-sx0)*(x-sx0) + (y-sx1)*(y-sx1));
			hij += 1.0/ rr * gausswt[gindex];
			gij += -log(sqrt(rr)) *gausswt[gindex];
		}
		d = (x-sx0)*unx + (y-sx1)*uny;
		hij *= -1.0/(4*pi)*gama_j*d;
		gij *= 1.0/(4*pi)*gama_j;

	//	std::cout << hij << " ++++ === +==+ ==++ ==+ =+ " << gij<<std::endl;

		ret += -hij*solution.u[j] + gij*solution.q[j];
	}
// 	std::cout << ret <<std::endl;

	return ret;
}

// calculate the h0 norm in the field.
// TODO: want to find an exact solution.
double Laplace2d::h0_norm()
{
	double gauss_num = 4;
	double gauss_pt[4] = { 0.86113631, -0.86113631, 0.33998014, -0.33998014 };
	double gauss_wt[4] = { 0.34785485, 0.34785485, 0.65214515, 0.65214515};

	double h0 = 0;
	for (int i = 0; i < bmesh.number_of_volume_elements(); ++ i)
	{
		const Mesh::Volume_element& ve = bmesh.find_volume_element(i);
		for (int g_id = 0; g_id < gauss_num; ++g_id )
		{
			h0 += 0.1; // TODO: 
		}
	}
	return h0;
}

void GAUSSELIMINATE(double **E,int M,double *RHS)
{
	using namespace std;
	int i;
	int k;
	for(k=0;k<M-1;k++)
	{
		double bmax=0.0;
		int ik;
		for(i=k;i<M;i++)
		{
			if(bmax<fabs(E[i][k]))
			{
				bmax=fabs(E[i][k]);
				ik=i;
			}
		}
		if(bmax<1.0e-10) 
		{
			cout<<"主元太小"<<endl;
			//system("pause"); // 主元太小
		}
		if(ik!=k)
		{
			double t;
			for(i=k;i<M;i++)
			{
				t=E[ik][i];
				E[ik][i]=E[k][i];
				E[k][i]=t;
			}
			t=RHS[ik];
			RHS[ik]=RHS[k];
			RHS[k]=t;
		}
		// 消元
		for(i=k+1;i<M;i++)
		{
			if(E[i][k]!=0)
			{
				double lk=E[i][k]/E[k][k];
				int j;
				for(j=k+1;j<M;j++)
				{
					E[i][j]=E[i][j]-lk*E[k][j];
				}
				RHS[i]=RHS[i]-lk*RHS[k];
			}
		}
	}
	if(fabs(E[M-1][M-1])<1.0e-10)
	{
		cout<<"主元太小"<<endl;
		//system("pause"); // 主元太小
	}

	// 消元法结束后开始回代
	RHS[M-1]=RHS[M-1]/E[M-1][M-1];
	for(i=M-2;i>=0;i--)
	{
		double s=0.0;
		int j;
		for(j=i+1;j<M;j++)
		{
			s=s+E[i][j]*RHS[j];
		}
		RHS[i]=(RHS[i]-s)/E[i][i];
	}
}

