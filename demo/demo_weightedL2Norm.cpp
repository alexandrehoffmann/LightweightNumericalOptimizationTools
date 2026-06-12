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
	
	using CG       = LNOT::ConjugateGradient<double>;
	using BisectLS = LNOT::BisectionLineSearch<double>;
	
	std::FILE* l2Out = std::fopen("newton_with_l2_norm.log", "w");
	
	LNOT::NewtonSolver<CG, BisectLS> newtonSolver1;
	newtonSolver1.setTol(1.0e-5);
	newtonSolver1.setOutput(l2Out);
	newtonSolver1.solve(func, grad, hessOp, precOp, N, x);
	
	fmt::print("Newton with L2 convergence criterion found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), newtonSolver1.getIterations(), newtonSolver1.getError(), newtonSolver1.getValue());
	
	// weighted L2 norm needs a scope, which contains the weights used by the norm.
	// When the scope is destroyed the data are freed
	using Criterion      = LNOT::WeightedL2NormConvergenceCriterion;
	using CriterionScope = typename Criterion::Scope<double, unsigned int>;
	
	std::FILE* weightedL2Out = std::fopen("newton_with_weighted_l2_norm.log", "w");
	
	CriterionScope scope(N);
	
	scope.getActiveWeights().front() = 1.0e8; // set the first weight to 10^{8}
	scope.getActiveWeights().back() = 1.0e-8; // set the last weight to 10^{-8}
	
	LNOT::NewtonSolver<CG, BisectLS, Criterion> newtonSolver2;
	newtonSolver2.setTol(1.0e-5);
	newtonSolver2.setOutput(weightedL2Out);
	newtonSolver2.solve(func, grad, hessOp, precOp, N, x);
	
	fmt::print("Newton with weighted L2 convergence criterion found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), newtonSolver2.getIterations(), newtonSolver2.getError(), newtonSolver2.getValue());
	
	std::fclose(l2Out);
	std::fclose(weightedL2Out);
  
	return EXIT_SUCCESS;
}
