#include <LNOT/Core.hpp>

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ranges.h>

#include <span>

int main()
{
	constexpr unsigned int N = 10;
	
 	const double A[N*N] = {
		0.42261004, 0.        , 0.        , 0.        , 0.        ,
		0.        , 0.        , 0.        , 0.        , 0.        ,
		0.        , 0.87694182, 0.        , 0.        , 0.        ,
		0.        , 0.        , 0.        , 0.        , 0.        ,
		0.        , 0.        , 0.4829434 , 0.        , 0.        ,
		0.        , 0.        , 0.        , 0.        , 0.        ,
		0.        , 0.        , 0.        , 0.43529297, 0.        ,
		0.        , 0.        , 0.        , 0.        , 0.        ,
		0.        , 0.        , 0.        , 0.        , 0.76719106,
		0.        , 0.        , 0.        , 0.        , 0.        ,
		0.        , 0.        , 0.        , 0.        , 0.        ,
		0.3279194 , 0.        , 0.        , 0.        , 0.        ,
		0.        , 0.        , 0.        , 0.        , 0.        ,
		0.        , 0.15173135, 0.        , 0.        , 0.        ,
		0.        , 0.        , 0.        , 0.        , 0.        ,
		0.        , 0.        , 0.10138197, 0.        , 0.        ,
		0.        , 0.        , 0.        , 0.        , 0.        ,
		0.        , 0.        , 0.        , 0.10743903, 0.        ,
		0.        , 0.        , 0.        , 0.        , 0.        ,
		0.        , 0.        , 0.        , 0.        , 0.00348401};
	const double b[N] = {
		0.92182983, 0.44848884, 0.67730341, 0.52299423, 0.31012989,
		0.03227012, 0.26467142, 0.43539253, 0.83558561, 0.16021716};

	double x[N];
	std::span<const double> x_view(x, N);
	
	LNOT::SymmetricDenseMatrixOp<double> Aop(A, N);
	LNOT::DiagonalPreconditionerOp<double> Bop(A, N);
	
	LNOT::ConjugateGradient<double> cg;
	cg.setOutput(stdout);
	cg.solve(Aop, Bop, b, BIC::fixed<unsigned int, N>, x);
	
	fmt::print("CG found : {:.2f} in {} iterations with a final error of {}\n", fmt::join(x_view, " "), cg.getIterations(), cg.getError());
	
	LNOT::LanczosSolver<double> lanczos;
	lanczos.setOutput(stdout);
	lanczos.solve(Aop, Bop, b, BIC::fixed<unsigned int, N>, x);
	
	fmt::print("Lanczos found : {:.2f} in {} iterations with a final error of {}\n", fmt::join(x_view, " "), lanczos.getIterations(), lanczos.getError());
	
	return EXIT_SUCCESS;
}
