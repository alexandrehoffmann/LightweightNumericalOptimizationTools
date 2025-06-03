#include <LNOT/Core.hpp>

#include <fmt/core.h>
#include <fmt/format.h>

#include <span>

int main()
{
	constexpr unsigned int N = 10;
	
 	const double A[N*N] = {4.06221475, 3.05878411, 3.28435976, 2.38055909, 2.41311328,
       2.7726784 , 3.69140663, 3.39073922, 3.76682296, 1.71341658,
       3.05878411, 4.06142791, 2.72437461, 2.37025037, 2.85672027,
       2.91818372, 3.71723443, 3.82002084, 3.73158745, 2.39528936,
       3.28435976, 2.72437461, 3.27629245, 2.17287238, 2.48131248,
       2.46303549, 2.87722868, 2.85118017, 3.45855898, 1.66731065,
       2.38055909, 2.37025037, 2.17287238, 2.20263528, 1.71011861,
       1.75354288, 2.88562711, 2.02965453, 2.49643072, 1.73761211,
       2.41311328, 2.85672027, 2.48131248, 1.71011861, 2.92850375,
       2.49253318, 2.82469135, 2.90761914, 3.08750857, 1.69948825,
       2.7726784 , 2.91818372, 2.46303549, 1.75354288, 2.49253318,
       2.99380146, 2.89355228, 3.03489047, 3.47917545, 1.68675963,
       3.69140663, 3.71723443, 2.87722868, 2.88562711, 2.82469135,
       2.89355228, 4.87088542, 3.45977699, 4.16030522, 2.49628993,
       3.39073922, 3.82002084, 2.85118017, 2.02965453, 2.90761914,
       3.03489047, 3.45977699, 4.08773201, 3.63006462, 1.87332754,
       3.76682296, 3.73158745, 3.45855898, 2.49643072, 3.08750857,
       3.47917545, 4.16030522, 3.63006462, 5.09014188, 2.2065241 ,
       1.71341658, 2.39528936, 1.66731065, 1.73761211, 1.69948825,
       1.68675963, 2.49628993, 1.87332754, 2.2065241 , 2.37528137};
	const double b[N] = {0.92182983, 0.44848884, 0.67730341, 0.52299423, 0.31012989,
       0.03227012, 0.26467142, 0.43539253, 0.83558561, 0.16021716};
  
  const double delta = 1;
  
	double x[N];
	std::span<const double> x_view(x, N);
	
	LNOT::SymmetricDenseMatrixOp<double> Aop(A, N);
	LNOT::LanczosTRSSolver<double> lanczosTrs;
	lanczosTrs.setOutput(stdout);
	lanczosTrs.solve(Aop, b, N, delta, x);
	
	fmt::print("Lanczos Trust Region found : {:.2f} in {} iterations with a final error of {:10.2e}\n", fmt::join(x_view, " "), lanczosTrs.getIterations(), lanczosTrs.getError());
	
	double r[N];
	std::span<const double> r_view(r, N);
	
	Aop(x, r);
	for (unsigned int i=0; i!=N; ++i) { r[i] += lanczosTrs.getLambda()*x[i] + b[i];  }
	
	fmt::print("KKT : \n");
	fmt::print("  * (A + lambda I)x + b = {:10.2e}\n", fmt::join(r_view, " "));
	fmt::print("  * |x| <= Delta ? {}\n", LNOT::BasicLinalg::norm(x, N) <= delta + lanczosTrs.getTolTR());
	fmt::print("  * lambda >= 0 ? {}\n", lanczosTrs.getLambda() >= 0);
	fmt::print("  * lambda*(|x| - Delta) = {:10.2e}\n", lanczosTrs.getLambda()*(LNOT::BasicLinalg::norm(x, N) - delta));
	Aop(x, r); 
	fmt::print("Predicted model reduction : {}\n", lanczosTrs.getModelReduction());
	fmt::print("Actual model reduction    : {}\n", 0.5*LNOT::BasicLinalg::inner(x, r, N) + LNOT::BasicLinalg::inner(x, b, N));
	
	LNOT::CoupledLanczosTRSSolver<double> coupledLanczosTr;
	coupledLanczosTr.setOutput(stdout);
	coupledLanczosTr.solve(Aop, b, N, delta, x);
	
	fmt::print("Lanczos Trust Region found : {:.2f} in {} iterations with a final error of {:10.2e}\n", fmt::join(x_view, " "), coupledLanczosTr.getIterations(), coupledLanczosTr.getError());
	
	Aop(x, r);
	for (unsigned int i=0; i!=N; ++i) { r[i] += coupledLanczosTr.getLambda()*x[i] + b[i];  }
	
	fmt::print("KKT : \n");
	fmt::print("  * (A + lambda I)x + b = {:10.2e}\n", fmt::join(r_view, " "));
	fmt::print("  * |x| <= Delta ? {}\n", LNOT::BasicLinalg::norm(x, N) <= delta + lanczosTrs.getTolTR());
	fmt::print("  * lambda >= 0 ? {}\n", coupledLanczosTr.getLambda() >= 0);
	fmt::print("  * lambda*(|x| - Delta) = {:10.2e}\n", coupledLanczosTr.getLambda()*(LNOT::BasicLinalg::norm(x, N) - delta));
	Aop(x, r); 
	fmt::print("Predicted model reduction : {}\n", coupledLanczosTr.getModelReduction());
	fmt::print("Actual model reduction    : {}\n", 0.5*LNOT::BasicLinalg::inner(x, r, N) + LNOT::BasicLinalg::inner(x, b, N));
	
	LNOT::TruncatedConjugateGradient<double> tcg;
	tcg.setOutput(stdout);
	tcg.solve(Aop, b, N, delta, x);
	
	fmt::print("Truncated CG found : {:.2f} in {} iterations with a final error of {:10.2e}\n", fmt::join(x_view, " "), tcg.getIterations(), tcg.getError());
	fmt::print("|x| <= Delta ? {}\n", LNOT::BasicLinalg::norm(x, N) <= delta + tcg.getTolTR());
	Aop(x, r); 
	fmt::print("Predicted model reduction : {}\n", tcg.getModelReduction());
	fmt::print("Actual model reduction    : {}\n", 0.5*LNOT::BasicLinalg::inner(x, r, N) + LNOT::BasicLinalg::inner(x, b, N));
	
	return EXIT_SUCCESS;
}
