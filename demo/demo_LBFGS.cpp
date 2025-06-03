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
	
	LNOT::TruncatedConjugateGradient<double> tcg;
	LNOT::LanczosTRSSolver<double> lanczosTrs;
	LNOT::CoupledLanczosTRSSolver<double> coupledLanczosTrs;
	
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

	auto sr1TR1 = LNOT::makeSR1Solver(tcg);
	sr1TR1.setTol(1.0e-11);
	sr1TR1.setOutput(stdout);
	sr1TR1.solve(func, grad, N, x);
	
	fmt::print("TR-SR1 found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), sr1TR1.getIterations(), sr1TR1.getError(), sr1TR1.getValue());

	auto sr1TR2 = LNOT::makeSR1Solver(lanczosTrs);
	sr1TR2.setTol(1.0e-11);
	sr1TR2.setOutput(stdout);
	sr1TR2.solve(func, grad, N, x);
	
	fmt::print("TR-SR1 found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), sr1TR2.getIterations(), sr1TR2.getError(), sr1TR2.getValue());

	auto sr1TR3 = LNOT::makeSR1Solver(coupledLanczosTrs);
	sr1TR3.setTol(1.0e-11);
	sr1TR3.setOutput(stdout);
	sr1TR3.solve(func, grad, N, x);
	
	fmt::print("TR-SR1 found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), sr1TR3.getIterations(), sr1TR3.getError(), sr1TR3.getValue());

	auto lsr1TR1 = LNOT::makeLSR1Solver(tcg, 5);
	lsr1TR1.setTol(1.0e-11);
	lsr1TR1.setOutput(stdout);
	lsr1TR1.solve(func, grad, N, x);
	
	fmt::print("TR-L-SR1 found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), lsr1TR1.getIterations(), lsr1TR1.getError(), lsr1TR1.getValue());

	auto lsr1TR2 = LNOT::makeLSR1Solver(lanczosTrs, 5);
	lsr1TR2.setTol(1.0e-11);
	lsr1TR2.setOutput(stdout);
	lsr1TR2.solve(func, grad, N, x);
	
	fmt::print("TR-L-SR1 found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), lsr1TR2.getIterations(), lsr1TR2.getError(), lsr1TR2.getValue());

	auto lsr1TR3 = LNOT::makeLSR1Solver(coupledLanczosTrs, 5);
	lsr1TR3.setTol(1.0e-11);
	lsr1TR3.setOutput(stdout);
	lsr1TR3.solve(func, grad, N, x);
	
	fmt::print("TR-L-SR1 found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), lsr1TR3.getIterations(), lsr1TR3.getError(), lsr1TR3.getValue());

	return EXIT_SUCCESS;
}
