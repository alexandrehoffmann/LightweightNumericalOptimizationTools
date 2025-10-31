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
	
	double diagA[N];
	
	LNOT::BasicLinalg::stridedCopy(A, N+1, N, diagA, 1u);
	
	const double delta = 1;
  
	double x[N];
	std::span<const double> x_view(x, N);
	
	LNOT::SymmetricDenseMatrixOp<double> Aop(A, N);
	LNOT::DiagonalPreconditionerOp<double> Bop(A, N);
	
	LNOT::LanczosTRSSolver<double> lanczosTrs;
	lanczosTrs.setOutput(stdout);
	lanczosTrs.solve(Aop, Bop, b, BIC::fixed<unsigned int, N>, delta, x);
	
	fmt::print("Lanczos Trust Region found : {:.2f} in {} iterations with a final error of {:10.2e}\n", fmt::join(x_view, " "), lanczosTrs.getIterations(), lanczosTrs.getError());
	
	double r[N];
	std::span<const double> r_view(r, N);
	
	Aop(x, r);
	for (unsigned int i=0; i!=N; ++i) { r[i] += lanczosTrs.getLambda()*A[i*N + i]*x[i] + b[i];  }
	
	fmt::print("KKT : \n");
	fmt::print("  * (A + lambda B)x + b = {:10.2e}\n", fmt::join(r_view, " "));
	fmt::print("  * |x|_B <= Delta ? {}\n", LNOT::BasicLinalg::weightedNorm(x, diagA, N) <= delta + lanczosTrs.getTolTR());
	fmt::print("  * lambda >= 0 ? {}\n", lanczosTrs.getLambda() >= 0);
	fmt::print("  * lambda*(|x| - Delta) = {:10.2e}\n", lanczosTrs.getLambda()*(LNOT::BasicLinalg::weightedNorm(x, diagA, N) - delta));
	Aop(x, r); 
	fmt::print("Predicted model reduction : {}\n", lanczosTrs.getModelReduction());
	fmt::print("Actual model reduction    : {}\n", 0.5*LNOT::BasicLinalg::inner(x, r, N) + LNOT::BasicLinalg::inner(x, b, N));
	
	LNOT::TruncatedConjugateGradient<double> tcg;
	tcg.setOutput(stdout);
	tcg.solve(Aop, Bop, b, BIC::fixed<unsigned int, N>, delta, x);
	
	fmt::print("Truncated CG found : {:.2f} in {} iterations with a final error of {:10.2e}\n", fmt::join(x_view, " "), tcg.getIterations(), tcg.getError());
	fmt::print("|x| <= Delta ? {}\n", LNOT::BasicLinalg::norm(x, N) <= delta + tcg.getTolTR());
	Aop(x, r); 
	fmt::print("Predicted model reduction : {}\n", tcg.getModelReduction());
	fmt::print("Actual model reduction    : {}\n", 0.5*LNOT::BasicLinalg::inner(x, r, N) + LNOT::BasicLinalg::inner(x, b, N));
	
	return EXIT_SUCCESS;
}
