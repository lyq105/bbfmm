#include "lapace2d.h"
#include "timer.h"
#include <string>
void constant_bem(const std::string& fname)
{
	std::string tec_fname = fname + ".plt";
	Timer tm;
	tm.start();
	Laplace2d laplace(fname);
	tm.stop();
	tm.diplayelapsedtime();
	laplace.initial_bc();
	tm.start();
	laplace.solve();
//	laplace.solve_lapack();
	tm.stop();	
	tm.diplayelapsedtime();
	laplace.print_solution( tec_fname );
}

int main(int argc, const char *argv[])
{
	constant_bem(std::string(argv[1]));
	return 0;
}
