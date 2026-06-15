#ifndef LNOT_NON_LINEAR_CONJUGATE_GRADIENT_IMPL_HPP
#define LNOT_NON_LINEAR_CONJUGATE_GRADIENT_IMPL_HPP

#include <LNOT/GeneralSolvers/NonLinearConjugateGradient.hpp>
#include <LNOT/LineSearches/BisectionLineSearch.hpp>
#include <LNOT/LineSearches/BacktrackingLineSearch.hpp>
#include <LNOT/LineSearches/NoLineSearch.hpp>

namespace LNOT
{
	
//// explicit template instanciations ////

//// for floats
extern template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::HESTENES_STIEFEL>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::FLETCHER_REEVES>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::POLAK_RIBIERE>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::FLETCHER>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::LIU_STOREY>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::DAI_YUAN>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::HAGER_ZHANG>;

extern template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::HESTENES_STIEFEL>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::FLETCHER_REEVES>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::POLAK_RIBIERE>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::FLETCHER>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::LIU_STOREY>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::DAI_YUAN>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::HAGER_ZHANG>;

extern template class NonLinearConjugateGradient<NoLineSearch<float>, NLCGUpdateStrategy::HESTENES_STIEFEL>;
extern template class NonLinearConjugateGradient<NoLineSearch<float>, NLCGUpdateStrategy::FLETCHER_REEVES>;
extern template class NonLinearConjugateGradient<NoLineSearch<float>, NLCGUpdateStrategy::POLAK_RIBIERE>;
extern template class NonLinearConjugateGradient<NoLineSearch<float>, NLCGUpdateStrategy::FLETCHER>;
extern template class NonLinearConjugateGradient<NoLineSearch<float>, NLCGUpdateStrategy::LIU_STOREY>;
extern template class NonLinearConjugateGradient<NoLineSearch<float>, NLCGUpdateStrategy::DAI_YUAN>;
extern template class NonLinearConjugateGradient<NoLineSearch<float>, NLCGUpdateStrategy::HAGER_ZHANG>;

//// for doubles
extern template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::HESTENES_STIEFEL>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::FLETCHER_REEVES>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::POLAK_RIBIERE>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::FLETCHER>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::LIU_STOREY>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::DAI_YUAN>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::HAGER_ZHANG>;

extern template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::HESTENES_STIEFEL>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::FLETCHER_REEVES>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::POLAK_RIBIERE>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::FLETCHER>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::LIU_STOREY>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::DAI_YUAN>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::HAGER_ZHANG>;

extern template class NonLinearConjugateGradient<NoLineSearch<double>, NLCGUpdateStrategy::HESTENES_STIEFEL>;
extern template class NonLinearConjugateGradient<NoLineSearch<double>, NLCGUpdateStrategy::FLETCHER_REEVES>;
extern template class NonLinearConjugateGradient<NoLineSearch<double>, NLCGUpdateStrategy::POLAK_RIBIERE>;
extern template class NonLinearConjugateGradient<NoLineSearch<double>, NLCGUpdateStrategy::FLETCHER>;
extern template class NonLinearConjugateGradient<NoLineSearch<double>, NLCGUpdateStrategy::LIU_STOREY>;
extern template class NonLinearConjugateGradient<NoLineSearch<double>, NLCGUpdateStrategy::DAI_YUAN>;
extern template class NonLinearConjugateGradient<NoLineSearch<double>, NLCGUpdateStrategy::HAGER_ZHANG>;

//// for long doubles
extern template class NonLinearConjugateGradient<BisectionLineSearch<long double>, NLCGUpdateStrategy::HESTENES_STIEFEL>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<long double>, NLCGUpdateStrategy::FLETCHER_REEVES>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<long double>, NLCGUpdateStrategy::POLAK_RIBIERE>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<long double>, NLCGUpdateStrategy::FLETCHER>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<long double>, NLCGUpdateStrategy::LIU_STOREY>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<long double>, NLCGUpdateStrategy::DAI_YUAN>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<long double>, NLCGUpdateStrategy::HAGER_ZHANG>;

extern template class NonLinearConjugateGradient<BacktrackingLineSearch<long double>, NLCGUpdateStrategy::HESTENES_STIEFEL>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<long double>, NLCGUpdateStrategy::FLETCHER_REEVES>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<long double>, NLCGUpdateStrategy::POLAK_RIBIERE>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<long double>, NLCGUpdateStrategy::FLETCHER>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<long double>, NLCGUpdateStrategy::LIU_STOREY>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<long double>, NLCGUpdateStrategy::DAI_YUAN>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<long double>, NLCGUpdateStrategy::HAGER_ZHANG>;

extern template class NonLinearConjugateGradient<NoLineSearch<long double>, NLCGUpdateStrategy::HESTENES_STIEFEL>;
extern template class NonLinearConjugateGradient<NoLineSearch<long double>, NLCGUpdateStrategy::FLETCHER_REEVES>;
extern template class NonLinearConjugateGradient<NoLineSearch<long double>, NLCGUpdateStrategy::POLAK_RIBIERE>;
extern template class NonLinearConjugateGradient<NoLineSearch<long double>, NLCGUpdateStrategy::FLETCHER>;
extern template class NonLinearConjugateGradient<NoLineSearch<long double>, NLCGUpdateStrategy::LIU_STOREY>;
extern template class NonLinearConjugateGradient<NoLineSearch<long double>, NLCGUpdateStrategy::DAI_YUAN>;
extern template class NonLinearConjugateGradient<NoLineSearch<long double>, NLCGUpdateStrategy::HAGER_ZHANG>;

//// method implementations ////

template<typename LineSearch, NLCGUpdateStrategy UpdateStrategy, class ConvergenceCriterion>  template<CFirstOrderOracle Oracle, typename ABool>
void NonLinearConjugateGradient<LineSearch, UpdateStrategy, ConvergenceCriterion>::solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value)
{
	using Oracle_Size = typename Oracle::Size;
	
	const Oracle_Size size = oracle.getNDims();
	
	if (not m_dk or m_workCapacity < size)
	{
		m_workCapacity = size;

		m_dk   = std::make_unique<Scalar[]>(m_workCapacity);
		m_yk   = std::make_unique<Scalar[]>(m_workCapacity);
		m_gk   = std::make_unique<Scalar[]>(m_workCapacity);
		m_gkp1 = std::make_unique<Scalar[]>(m_workCapacity);
	}
	if (not solveInPlace) { std::fill(x, x + size, 0); }
	
	m_innerIts.clear();
	
	constexpr FPComparator<Scalar> cmp;
	
	oracle.setCurrentPoint(x);
	oracle.getGradient(m_gk.get());
	
	m_criterion.init(m_gk.get(), size);
	
	#pragma omp simd
	for (Size i=0; i!=size; ++i) { m_dk[i] = -m_gk[i]; }
	
	m_fx = oracle.getValue();
	m_residual = m_criterion.getInitialResidual();
	
	const Scalar relTol   = m_criterion.getRelTol(m_relTol);
	const Scalar absTol   = m_criterion.getAbsTol(m_absTol);
	
	if (m_out != nullptr) { fmt::println(m_out, "#Non Linear CG method\n#Iteration f(x) residual relative_tol absolute_tol"); }
	
	m_info = Info::FAILURE;
	for (m_nIt=0;m_nIt!=m_maxIt; ++m_nIt)
	{
		if (m_out) { fmt::println(m_out, "{} {:10.2e} {:10.2e} {:10.2e} {:10.2e}", m_nIt, m_fx, m_residual, relTol, absTol); std::fflush(m_out); }
		if (m_residual < relTol or m_residual < absTol) { m_info = Info::SUCCESS; break; }
		
		const Scalar alpha = m_lineSearch.solve(x, m_fx, m_gk.get(), m_dk.get(), oracle);
		
		if (not cmp.isDefPositive(alpha*BasicLinalg::norm(m_dk.get(), size))) { m_info = Info::BREAKDOWN; break; }
		
		m_innerIts.push_back(1);
		
		BasicLinalg::axpy(alpha, m_dk.get(), size, x);
		
		oracle.setCurrentPoint(x);
		oracle.getGradient(m_gkp1.get());
		
		m_fx = oracle.getValue();
		m_residual = m_criterion.getResidual(m_gk.get(), size);
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_yk[i] = m_gkp1[i] - m_gk[i]; }
		
		const Scalar beta = getBeta(size);
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_dk[i] = -m_gkp1[i] + beta*m_dk[i]; }
		
		std::copy(m_gkp1.get(), m_gkp1.get() + size, m_gk.get());
	}
}

template<typename LineSearch, NLCGUpdateStrategy UpdateStrategy, class ConvergenceCriterion> 
auto NonLinearConjugateGradient<LineSearch, UpdateStrategy, ConvergenceCriterion>::getBeta(const Size size) -> Scalar
{
	if      constexpr (UpdateStrategy == NLCGUpdateStrategy::HESTENES_STIEFEL) { return BasicLinalg::inner(m_gkp1.get(), m_yk.get(), size) /  BasicLinalg::inner(m_dk.get(), m_yk.get(), size); }
	else if constexpr (UpdateStrategy == NLCGUpdateStrategy::FLETCHER_REEVES)  { return BasicLinalg::squaredNorm(m_gkp1.get(), size)       /  BasicLinalg::squaredNorm(m_gk.get(), size); }
	else if constexpr (UpdateStrategy == NLCGUpdateStrategy::POLAK_RIBIERE)    { return BasicLinalg::inner(m_gkp1.get(), m_yk.get(), size) /  BasicLinalg::squaredNorm(m_gk.get(), size); }
	else if constexpr (UpdateStrategy == NLCGUpdateStrategy::FLETCHER)         { return BasicLinalg::squaredNorm(m_gkp1.get(), size)       / -BasicLinalg::inner(m_dk.get(), m_gk.get(), size); }
	else if constexpr (UpdateStrategy == NLCGUpdateStrategy::LIU_STOREY)       { return BasicLinalg::inner(m_gkp1.get(), m_yk.get(), size) / -BasicLinalg::inner(m_dk.get(), m_gk.get(), size); }
	else if constexpr (UpdateStrategy == NLCGUpdateStrategy::DAI_YUAN)         { return BasicLinalg::squaredNorm(m_gkp1.get(), size)       /  BasicLinalg::inner(m_dk.get(), m_yk.get(), size); }
	else                                                                       
	{ 
		const Scalar sqNormY = BasicLinalg::squaredNorm(m_yk.get(), size);
		const Scalar yDotG   = BasicLinalg::inner(m_dk.get(), m_gk.get(), size);
		const Scalar beta1   = BasicLinalg::inner(m_yk.get(), m_gkp1.get(), size) / yDotG;
		const Scalar beta2   = -2*(sqNormY / (yDotG*yDotG))*BasicLinalg::inner(m_dk.get(), m_gkp1.get(), size);
		return beta1 + beta2;
	}
}

} // namespace LNOT

#endif // LNOT_NON_LINEAR_CONJUGATE_GRADIENT_IMPL_HPP
