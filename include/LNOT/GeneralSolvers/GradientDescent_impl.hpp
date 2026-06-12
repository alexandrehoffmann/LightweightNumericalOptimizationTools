#ifndef LNOT_GRADIENT_DESCENT_IMPL_HPP
#define LNOT_GRADIENT_DESCENT_IMPL_HPP

#include <LNOT/GeneralSolvers/GradientDescent.hpp>
#include <LNOT/LineSearches/BisectionLineSearch.hpp>
#include <LNOT/LineSearches/BacktrackingLineSearch.hpp>
#include <LNOT/LineSearches/NoLineSearch.hpp>

namespace LNOT
{

//// explicit template instanciations ////

extern template class GradientDescent< BisectionLineSearch<float> >;
extern template class GradientDescent< BacktrackingLineSearch<float> >;
extern template class GradientDescent< NoLineSearch<float> >;

extern template class GradientDescent< BisectionLineSearch<double> >;
extern template class GradientDescent< BacktrackingLineSearch<double> >;
extern template class GradientDescent< NoLineSearch<double> >;

extern template class GradientDescent< BisectionLineSearch<long double> >;
extern template class GradientDescent< BacktrackingLineSearch<long double> >;
extern template class GradientDescent< NoLineSearch<long double> >;

//// method implementations ////

template<typename LineSearch, typename ConvergenceCriterion> template<CFirstOrderOracle Oracle, typename ABool> 
void GradientDescent<LineSearch, ConvergenceCriterion>::solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value)
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
		
	m_innerIts.clear();
	
	constexpr ConvergenceCriterion criterion;
	constexpr FPComparator<Scalar> cmp;
	
	oracle.setCurrentPoint(x);
	oracle.getGradient(m_gk.get());
	
	m_fx       = oracle.getValue();
	m_residual = criterion.getResidual(m_gk.get(), size);
	
	const Scalar relTol = criterion.getRelTol(m_relTol, m_residual);
	const Scalar absTol = criterion.getAbsTol(m_absTol);
	
	if (m_out != nullptr) { fmt::println(m_out, "#Gradient descent method\n#Iteration f(x) residual relative_tol absolute_tol"); }
	
	m_info = Info::FAILURE;
	for (m_nIt=0;m_nIt!=m_maxIt; ++m_nIt)
	{
		if (m_out) { fmt::println(m_out, "{} {:10.2e} {:10.2e} {:10.2e} {:10.2e}", m_nIt, m_fx, m_residual, relTol, absTol); std::fflush(m_out); }
		if (m_residual < relTol or m_residual < absTol) { m_info = Info::SUCCESS; break; }
		
		#pragma omp simd
		for (Size j=0; j!=size; ++j) { m_sk[j] = -m_gk[j]; }
		m_innerIts.push_back(1);
		
		////
		const Scalar step_len = m_lineSearch.solve(x, m_fx, m_gk.get(), m_sk.get(), oracle);
		
		if (not cmp.isDefPositive(step_len*BasicLinalg::norm(m_sk.get(), size))) { m_info = Info::BREAKDOWN; break; }
		
		BasicLinalg::axpy(step_len, m_sk.get(), size, x);
		
		oracle.setCurrentPoint(x);
		oracle.getGradient(m_gk.get());
		
		m_fx = oracle.getValue();
		m_residual = criterion.getResidual(m_gk.get(), size);
	}
}

} // namespace LNOT

#endif // LNOT_GRADIENT_DESCENT_IMPL_HPP
