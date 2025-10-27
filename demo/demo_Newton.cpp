#include <LNOT/Core.hpp>

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ranges.h>

#include <span>

int main()
{
	using Size = unsigned int;
	constexpr Size N = 10;
	
	const auto func = [](const double* x)
	{ 
		double y(0);
		for (Size i=0; i!=Size(N-1); ++i)
		{
			y += 100*(x[i+1] - x[i]*x[i])*(x[i+1] - x[i]*x[i]) + (1 - x[i])*(1 - x[i]);
		}
		return y;
	};
	const auto grad = [](const double* x, double* g)
	{
		g[0] = -400*x[0]*(x[1] - x[0]*x[0]) + 2*(x[0] - 1);
		for (Size i=1; i!=Size(N-1); ++i)
		{
			g[i] = 200*(x[i] - x[i-1]*x[i-1]) - 400*x[i]*(x[i+1] - x[i]*x[i]) + 2*(x[i] - 1);
		}
		g[N-1] = 200*(x[N-1] - x[N-2]*x[N-2]);
	};
	const auto hessOp = [](const double* x, const double* d, double* Hd) 
	{
		std::fill(Hd, Hd + N, 0);
	
		Hd[0] += (2 - 400*(x[1] - 3*x[0]*x[0]))*d[0];
		Hd[1] += -400*x[0]*d[0];
		
		for (Size i=1; i!=Size(N-1); ++i)
		{
			Hd[i-1] += -400*x[i-1]*d[i]; 
			Hd[i]   += (202 - 400*(x[i+1] - 3*x[i]*x[i]))*d[i];
			Hd[i+1] += -400*x[i]*d[i];
		}
		
		Hd[N-2] += -400*x[N-2]*d[N-1];
		Hd[N-1] += 200*d[N-1];
	};
	const auto precOp = [](const double* x, const double* d, double* Hd) 
	{
		std::fill(Hd, Hd + N, 0);
	
		Hd[0] += d[0] / (2 - 400*(x[1] - 3*x[0]*x[0]));
		for (Size i=1; i!=Size(N-1); ++i)
		{
			Hd[i] += d[i] / (202 - 400*(x[i+1] - 3*x[i]*x[i]));
		}
		Hd[N-1] += d[N-1] / 200;
	};
	
	double x[N];
	std::span<const double> x_view(x, N);
	
	LNOT::ConjugateGradient<double> cg;
	LNOT::BisectionLineSearch<double> bisectLs;
	LNOT::BacktrackingLineSearch<double> backtrackLs;
	
	LNOT::TruncatedConjugateGradient<double> tcg;
	LNOT::LanczosTRSSolver<double> lanczosTrs;
	LNOT::CoupledLanczosTRSSolver<double> coupledLanczosTrs;
	
	std::FILE* newtonBisectLsOut       = std::fopen("newton_bisect_ls.log", "w");
	std::FILE* newtonBacktrackLsOut    = std::fopen("newton_backtrack_ls.log", "w");
	std::FILE* newtonTCGOut            = std::fopen("newton_tcg.log", "w");
	std::FILE* newtonLanczosOut        = std::fopen("newton_lanczos.log", "w");
	std::FILE* newtonCoupledLanczosOut = std::fopen("newton_coupled_lanczos.log", "w");
	
	auto newtonSolver1 = LNOT::makeNewtonSolver(cg, bisectLs);
	newtonSolver1.setOutput(newtonBisectLsOut);
	newtonSolver1.solve(func, grad, hessOp, precOp, N, x);
	
	fmt::print("Newton with Bisection LineSrearch found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), newtonSolver1.getIterations(), newtonSolver1.getError(), newtonSolver1.getValue());
	
	auto newtonSolver2 = LNOT::makeNewtonSolver(cg, backtrackLs);
	newtonSolver2.setOutput(newtonBacktrackLsOut);
	newtonSolver2.solve(func, grad, hessOp, precOp, N, x);
	
	fmt::print("Newton with Bisection Backtracking LineSearch found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), newtonSolver2.getIterations(), newtonSolver2.getError(), newtonSolver2.getValue());
	
	auto newtonSolver3 = LNOT::makeNewtonSolver(tcg);
	newtonSolver3.setOutput(newtonTCGOut);
	newtonSolver3.solve(func, grad, hessOp, N, x);
	
	fmt::print("Trut Region Newton with TCG found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), newtonSolver3.getIterations(), newtonSolver3.getError(), newtonSolver3.getValue());
	
	auto newtonSolver4 = LNOT::makeNewtonSolver(lanczosTrs);
	newtonSolver4.setOutput(newtonLanczosOut);
	newtonSolver4.solve(func, grad, hessOp, N, x);
	
	fmt::print("Trut Region Newton with Lanczos solver found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), newtonSolver4.getIterations(), newtonSolver4.getError(), newtonSolver4.getValue());
	
	auto newtonSolver5 = LNOT::makeNewtonSolver(coupledLanczosTrs);
	newtonSolver5.setOutput(newtonCoupledLanczosOut);
	newtonSolver5.solve(func, grad, hessOp, N, x);
	
	fmt::print("Trut Region Newton with coupled Lanczos solver found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), newtonSolver5.getIterations(), newtonSolver5.getError(), newtonSolver5.getValue());

	std::fclose(newtonBisectLsOut);
	std::fclose(newtonBacktrackLsOut);
	std::fclose(newtonTCGOut);
	std::fclose(newtonLanczosOut);
  std::fclose(newtonCoupledLanczosOut);
  
	return EXIT_SUCCESS;
}
