#ifndef LNOT_BISECTION_LINE_SEARCH_IMPL_HPP
#define LNOT_BISECTION_LINE_SEARCH_IMPL_HPP

#include <LNOT/LineSearches/BisectionLineSearch.hpp>
#include <LNOT/BasicLinalg/BasicLinalg.hpp>
#include <LNOT/FloatingPoint/FPComparator.hpp>

#include <fmt/core.h>
#include <fmt/format.h>

namespace LNOT
{
	
//// explicit template instanciations ////

extern template class BisectionLineSearch<float>;
extern template class BisectionLineSearch<double>;
extern template class BisectionLineSearch<long double>;

//// method implementations ////
	
template<typename T>
BisectionLineSearch<T>::BisectionLineSearch(const Size maxIt, const Scalar firstWolfConditionConst, const Scalar secondWolfConditionConst)
	: Base(maxIt)
	, m_firstWolfConditionConst(firstWolfConditionConst)
	, m_secondWolfConditionConst(secondWolfConditionConst)
{
}
	
template<typename T> template<CFirstOrderOracle Oracle>
auto BisectionLineSearch<T>::solveImpl(const Scalar* x, const Scalar& fx, const Scalar* gx, const Scalar* s, Oracle& oracle) -> Scalar
{
	using Oracle_Size = typename Oracle::Size;
	
	const Scalar inf = NumTraits<Scalar>::infinity;
	
	const Oracle_Size size     = oracle.getNDims();
	const Scalar      sDotGrad = BasicLinalg::inner(s, gx, size);
	
	constexpr FPComparator<Scalar> cmp;
	
	if (not m_newX or m_workCapacity < size)
	{
		m_workCapacity = size;
		m_newX    = std::make_unique<Scalar[]>(m_workCapacity);
		m_newGrad = std::make_unique<Scalar[]>(m_workCapacity);
	}
	
	const Scalar normS = BasicLinalg::norm(s, size);
	
	Scalar alpha    (1);
	Scalar alpha_min{};
	Scalar alpha_max(inf);
	
	m_info = Info::FAILURE;
	if (m_out != nullptr) { fmt::println(m_out, "#Bisection LineSearch : \n#Iteration alpha"); }
	for (m_nIt=0; m_nIt!=m_maxIt; ++m_nIt)
	{		
		if (m_out != nullptr) { fmt::println(m_out, "{} {:10.2e} {:10.2e} {:10.2e}", m_nIt, alpha, alpha_min, alpha_max); std::fflush(m_out); }
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_newX[i] = x[i] + alpha*s[i]; }
		
		oracle.setCurrentPoint(m_newX.get());
		oracle.getGradient(m_newGrad.get());
		
		if (oracle.isFeasible() and cmp.isApproxEq(alpha_max, alpha_min)) { break; } // (alpha_min, alpha_max) is empty
		
		const Scalar fx_new = oracle.getValue();
		const Scalar new_sDotGrad = BasicLinalg::inner(s, m_newGrad.get(), size);
		
		const bool isStepTooShort = cmp.isDefLessThan(new_sDotGrad, m_secondWolfConditionConst*sDotGrad);
		const bool isStepTooLong  = cmp.isDefGreaterThan(fx_new, fx + m_firstWolfConditionConst*alpha*sDotGrad);
				
		if (not oracle.isFeasible())
		{
			alpha_max = alpha;
			alpha = Scalar(0.5)*(alpha_max + alpha_min);
		}
		else if (isStepTooLong)
		{
			alpha_max = alpha;
			alpha = Scalar(0.5)*(alpha_max + alpha_min);
		}
		else if (isStepTooShort)
		{
			alpha_min = alpha;
			alpha = (alpha_max == inf) ? Scalar(2)*alpha : Scalar(0.5)*(alpha_max + alpha_min);
		}
		else
		{
			m_info = Info::SUCCESS;
			break;
		}
	}
	return alpha;
}
	
} // namespace LNOT

#endif // LNOT_BISECTION_LINE_SEARCH_IMPL_HPP

