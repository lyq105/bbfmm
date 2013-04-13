#include "lapace2d.h"
#include "timer.h"

void constant_bem()
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

int main(int argc, const char *argv[])
{
	constant_bem();
	return 0;
}
