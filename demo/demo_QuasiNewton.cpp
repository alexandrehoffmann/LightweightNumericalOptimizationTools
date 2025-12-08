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
	
	double x[N];
	std::span<const double> x_view(x, N);
	
	using BisectLs    = LNOT::BisectionLineSearch<double>;
	using BacktrackLs = LNOT::BacktrackingLineSearch<double>;

	using TCG        = LNOT::TruncatedConjugateGradient<double>;
	using LanczosTrs = LNOT::LanczosTRSSolver<double>;
	
	std::FILE* nlcgBisectLsOut       = std::fopen("nlcg_bisect_ls.log", "w");
	std::FILE* nlcgBacktrackLsOut    = std::fopen("nlcg_backtrack_ls.log", "w");
	std::FILE* bfgsBisectLsOut       = std::fopen("bfgs_bisect_ls.log", "w");
	std::FILE* bfgsBacktrackLsOut    = std::fopen("bfgs_backtrack_ls.log", "w");
	std::FILE* lbfgsBisectLsOut      = std::fopen("lbfgs_bisect_ls.log", "w");
	std::FILE* lbfgsBacktrackLsOut   = std::fopen("lbfgs_backtrack_ls.log", "w");
	std::FILE* sr1TCGOut             = std::fopen("sr1_tcg.log", "w");
	std::FILE* sr1LanczosOut         = std::fopen("sr1_lanczos.log", "w");
	std::FILE* lsr1TCGOut            = std::fopen("lsr1_tcg.log", "w");
	std::FILE* lsr1LanczosOut        = std::fopen("lsr1_lanczos.log", "w");
	
	LNOT::NonLinearConjugateGradient<BisectLs> nlcg1;
	nlcg1.setTol(1.0e-11);
	nlcg1.setOutput(nlcgBisectLsOut);
	nlcg1.solve(func, grad, N, x);
	
	fmt::print("NLCG found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), nlcg1.getIterations(), nlcg1.getError(), nlcg1.getValue());

	LNOT::NonLinearConjugateGradient<BacktrackLs> nlcg2;
	nlcg2.setTol(1.0e-11);
	nlcg2.setOutput(nlcgBacktrackLsOut);
	nlcg2.solve(func, grad, N, x);
	
	fmt::print("NLCG found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), nlcg2.getIterations(), nlcg2.getError(), nlcg2.getValue());
	
	LNOT::BFGS<BisectLs> bfgs1;
	bfgs1.setTol(1.0e-11);
	bfgs1.setOutput(bfgsBisectLsOut);
	bfgs1.solve(func, grad, N, x);
	
	fmt::print("BFGS found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), bfgs1.getIterations(), bfgs1.getError(), bfgs1.getValue());

	LNOT::BFGS<BacktrackLs> bfgs2;
	bfgs2.setTol(1.0e-11);
	bfgs2.setOutput(bfgsBacktrackLsOut);
	bfgs2.solve(func, grad, N, x);
	
	fmt::print("BFGS found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), bfgs2.getIterations(), bfgs2.getError(), bfgs2.getValue());
	
	LNOT::LBFGS<BisectLs> lbfgs1;
	lbfgs1.setTol(1.0e-11);
	lbfgs1.setOutput(lbfgsBisectLsOut);
	lbfgs1.solve(func, grad, N, x);
	
	fmt::print("L-BFGS found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), lbfgs1.getIterations(), lbfgs1.getError(), lbfgs1.getValue());

	LNOT::LBFGS<BacktrackLs> lbfgs2;
	lbfgs2.setTol(1.0e-11);
	lbfgs2.setOutput(lbfgsBacktrackLsOut);
	lbfgs2.solve(func, grad, N, x);
	
	fmt::print("L-BFGS found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), lbfgs2.getIterations(), lbfgs2.getError(), lbfgs2.getValue());
	
	LNOT::SR1TrustRegionSolver<TCG> sr1TR1;
	sr1TR1.setTol(1.0e-11);
	sr1TR1.setOutput(sr1TCGOut);
	sr1TR1.solve(func, grad, N, x);
	
	fmt::print("TR-SR1 found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), sr1TR1.getIterations(), sr1TR1.getError(), sr1TR1.getValue());
	
	LNOT::SR1TrustRegionSolver<LanczosTrs> sr1TR2;
	sr1TR2.setTol(1.0e-11);
	sr1TR2.setOutput(sr1LanczosOut);
	sr1TR2.solve(func, grad, N, x);
	
	fmt::print("TR-SR1 found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), sr1TR2.getIterations(), sr1TR2.getError(), sr1TR2.getValue());
	
	LNOT::LSR1TrustRegionSolver<TCG> lsr1TR1;
	lsr1TR1.setTol(1.0e-11);
	lsr1TR1.setOutput(lsr1TCGOut);
	lsr1TR1.solve(func, grad, N, x);
	
	fmt::print("TR-L-SR1 found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), lsr1TR1.getIterations(), lsr1TR1.getError(), lsr1TR1.getValue());
	
	LNOT::LSR1TrustRegionSolver<LanczosTrs> lsr1TR2;
	lsr1TR2.setTol(1.0e-11);
	lsr1TR2.setOutput(lsr1LanczosOut);
	lsr1TR2.solve(func, grad, N, x);
	
	fmt::print("TR-L-SR1 found : {:.2f} in {} iterations with a final error of {} and f(x) = {}\n", fmt::join(x_view, " "), lsr1TR2.getIterations(), lsr1TR2.getError(), lsr1TR2.getValue());

	std::fclose(bfgsBisectLsOut);
	std::fclose(bfgsBacktrackLsOut);
	std::fclose(lbfgsBisectLsOut);
	std::fclose(lbfgsBacktrackLsOut);
	std::fclose(sr1TCGOut);
	std::fclose(sr1LanczosOut);
	std::fclose(lsr1TCGOut);
	std::fclose(lsr1LanczosOut);
	
	return EXIT_SUCCESS;
}
