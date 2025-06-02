#include <LNOT/Core.hpp>

#include <fmt/core.h>
#include <fmt/format.h>

#include <span>

int main()
{
	using Size = unsigned int;
	constexpr Size N = 10;
	
	const auto func = [](const double* __restrict__ x)
	{ 
		double y(0);
		for (Size i=0; i!=Size(N-1); ++i)
		{
			y += 100*(x[i+1] - x[i]*x[i])*(x[i+1] - x[i]*x[i]) + (1 - x[i])*(1 - x[i]);
		}
		return y;
	};
	const auto grad = [](const double* __restrict__ x, double* __restrict__ g)
	{
		g[0] = -400*x[0]*(x[1] - x[0]*x[0]) + 2*(x[0] - 1);
		for (Size i=1; i!=Size(N-1); ++i)
		{
			g[i] = 200*(x[i] - x[i-1]*x[i-1]) - 400*x[i]*(x[i+1] - x[i]*x[i]) + 2*(x[i] - 1);
		}
		g[N-1] = 200*(x[N-1] - x[N-2]*x[N-2]);
	};
	const auto hessOp = [](const double* __restrict__ x, const double* __restrict__ d, double* __restrict__ Hd) 
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
	
	double x[N];
	std::span<const double> x_view(x, N);
	
	LNOT::ConjugateGradient<double> cg;
	LNOT::BisectionLineSearch<double> bisectLs;
	LNOT::BacktrackingLineSearch<double> backtrackLs;
	
	std::FILE* bisectLsOut = std::fopen("bisectLs.log", "w");
	std::FILE* backtrackLsOut = std::fopen("backtrackLs.log", "w");
	
	bisectLs.setOutput(bisectLsOut);
	backtrackLs.setOutput(backtrackLsOut);
	
	auto newtonSolver1 = LNOT::makeNewtonSolver(cg, bisectLs);
	newtonSolver1.setOutput(stdout);
	newtonSolver1.solve(func, grad, hessOp, N, x);
	
	fmt::print("Newton found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), newtonSolver1.getIterations(), newtonSolver1.getError(), newtonSolver1.getValue());
	for (Size it=0; it!=newtonSolver1.getIterations(); ++it)
	{
		fmt::print("  * Iteration {} computed in {} CG steps\n", it, newtonSolver1.getInnerIterations(it));
	}
	
	auto newtonSolver2 = LNOT::makeNewtonSolver(cg, backtrackLs);
	newtonSolver2.setOutput(stdout);
	newtonSolver2.solve(func, grad, hessOp, N, x);
	
	fmt::print("Newton found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), newtonSolver2.getIterations(), newtonSolver2.getError(), newtonSolver2.getValue());
	for (Size it=0; it!=newtonSolver2.getIterations(); ++it)
	{
		fmt::print("  * Iteration {} computed in {} CG steps\n", it, newtonSolver2.getInnerIterations(it));
	}
	
	std::fclose(bisectLsOut);
	std::fclose(backtrackLsOut);

	return EXIT_SUCCESS;
}
