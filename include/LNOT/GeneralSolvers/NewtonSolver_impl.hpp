#ifndef LNOT_NEWTON_SOLVER_IMPL_HPP
#define LNOT_NEWTON_SOLVER_IMPL_HPP

#include <LNOT/GeneralSolvers/NewtonSolver.hpp>
#include <LNOT/LineSearches/BisectionLineSearch.hpp>
#include <LNOT/LineSearches/BacktrackingLineSearch.hpp>
#include <LNOT/LineSearches/NoLineSearch.hpp>
#include <LNOT/LinearSolvers/ConjugateGradient.hpp>
#include <LNOT/LinearSolvers/LanczosSolver.hpp>

#include <fmt/core.h>
#include <fmt/format.h>

namespace LNOT
{

//// explicit template instanciations ////

// explicit instanciation for float
extern template class NewtonSolver<ConjugateGradient<float>, BisectionLineSearch<float> >;
extern template class NewtonSolver<LanczosSolver<float>,     BisectionLineSearch<float> >;

extern template class NewtonSolver<ConjugateGradient<float>, BacktrackingLineSearch<float> >;
extern template class NewtonSolver<LanczosSolver<float>,     BacktrackingLineSearch<float> >;

extern template class NewtonSolver<ConjugateGradient<float>, NoLineSearch<float> >;
extern template class NewtonSolver<LanczosSolver<float>,     NoLineSearch<float> >;

// explicit instanciation for double
extern template class NewtonSolver<ConjugateGradient<double>, BisectionLineSearch<double> >;
extern template class NewtonSolver<LanczosSolver<double>,     BisectionLineSearch<double> >;

extern template class NewtonSolver<ConjugateGradient<double>, BacktrackingLineSearch<double> >;
extern template class NewtonSolver<LanczosSolver<double>,     BacktrackingLineSearch<double> >;

extern template class NewtonSolver<ConjugateGradient<double>, NoLineSearch<double> >;
extern template class NewtonSolver<LanczosSolver<double>,     NoLineSearch<double> >;

// explicit instanciation for long double
extern template class NewtonSolver<ConjugateGradient<long double>, BisectionLineSearch<long double> >;
extern template class NewtonSolver<LanczosSolver<long double>,     BisectionLineSearch<long double> >;

extern template class NewtonSolver<ConjugateGradient<long double>, BacktrackingLineSearch<long double> >;
extern template class NewtonSolver<LanczosSolver<long double>,     BacktrackingLineSearch<long double> >;

extern template class NewtonSolver<ConjugateGradient<long double>, NoLineSearch<long double> >;
extern template class NewtonSolver<LanczosSolver<long double>,     NoLineSearch<long double> >;

//// method implementations ////

template<typename LinSolver, typename LineSearch, typename ConvergenceCriterion> template<CSecondOrderOracle Oracle, typename ABool> 
void NewtonSolver<LinSolver,LineSearch,ConvergenceCriterion>::solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value)
{
	using Oracle_Size = typename Oracle::Size;
	
	const Oracle_Size size = oracle.getNDims();
	
	if (not m_gk or m_workCapacity < size)
	{
		m_workCapacity = size;
		m_gk = std::make_unique<Scalar[]>(m_workCapacity);
		m_sk = std::make_unique<Scalar[]>(m_workCapacity);
	}
	if (not solveInPlace) { std::fill(x, x + size, 0); }
	
	const auto Hk    = [&oracle] (const Scalar* d, Scalar* Hd)    -> void { oracle.getHessianProd(d, Hd);  };
	const auto invBk = [&oracle] (const Scalar* d, Scalar* invBd) -> void { oracle.applyPrecond(d, invBd); };
	
	m_innerIts.clear();
	
	constexpr ConvergenceCriterion criterion;
	constexpr FPComparator<Scalar> cmp;
	
	oracle.setCurrentPoint(x);
	oracle.getGradient(m_gk.get());
	
	m_fx = oracle.getValue();
	m_residual = criterion.getResidual(m_gk.get(), size);
	
	const Scalar relTol = criterion.getRelTol(m_relTol, m_residual);
	const Scalar absTol = criterion.getAbsTol(m_absTol);
	
	if (m_out != nullptr) { fmt::println(m_out, "#Newton method\n#Iteration f(x) residual relative_tol absolute_tol"); }
	
	m_info = Info::FAILURE;
	for (m_nIt=0;m_nIt!=m_maxIt; ++m_nIt)
	{				
		if (m_out) { fmt::println(m_out, "{} {:10.2e} {:10.2e} {:10.2e} {:10.2e}", m_nIt, m_fx, m_residual, relTol, absTol); std::fflush(m_out); }
		if (m_residual < relTol or m_residual < absTol) { m_info = Info::SUCCESS; break; }
		
		m_linsSolver.solve(Hk, invBk, m_gk.get(), size, m_sk.get()); 
		
		if (m_linsSolver.getInfo() == LinSolver::Info::NEGATIVE_CURVATURE) 
		{ 
			#pragma omp simd
			for (Size i=0; i!=size; ++i) { m_sk[i] = -m_gk[i]; }
		}
		m_innerIts.push_back(m_linsSolver.getIterations());
		
		const Scalar alpha = m_lineSearch.solve(x, m_fx, m_gk.get(), m_sk.get(), oracle);
		
		if (not cmp.isDefPositive(alpha)) { m_info = Info::BREAKDOWN; break; }
		
		BasicLinalg::axpy(alpha, m_sk.get(), size, x);
		
		oracle.setCurrentPoint(x);
		oracle.getGradient(m_gk.get());
		
		m_fx = oracle.getValue();
		m_residual = criterion.getResidual(m_gk.get(), size);
	}
}

} // namespace LNOT

#endif // LNOT_NEWTON_SOLVER_IMPL_HPP
