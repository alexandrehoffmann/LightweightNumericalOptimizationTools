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
	
	double x[N];
	std::span<const double> x_view(x, N);
	
	LNOT::BisectionLineSearch<double> bisectLs;
	LNOT::BacktrackingLineSearch<double> backtrackLs;
	
	std::FILE* bisectLsOut = std::fopen("bisectLs.log", "w");
	std::FILE* backtrackLsOut = std::fopen("backtrackLs.log", "w");
	
	bisectLs.setOutput(bisectLsOut);
	backtrackLs.setOutput(backtrackLsOut);
	
	auto lbfgs1 = LNOT::makeLBFGS(bisectLs, 5);
	lbfgs1.setTol(1.0e-11);
	lbfgs1.setOutput(stdout);
	lbfgs1.solve(func, grad, N, x);
	
	fmt::print("LBFGS found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), lbfgs1.getIterations(), lbfgs1.getError(), lbfgs1.getValue());

	auto lbfgs2 = LNOT::makeLBFGS(backtrackLs, 5);
	lbfgs2.setTol(1.0e-11);
	lbfgs2.setOutput(stdout);
	lbfgs2.solve(func, grad, N, x);
	
	fmt::print("LBFGS found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), lbfgs2.getIterations(), lbfgs2.getError(), lbfgs2.getValue());
	
	std::fclose(bisectLsOut);
	std::fclose(backtrackLsOut);

	return EXIT_SUCCESS;
}
