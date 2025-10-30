#ifndef LNOT_BISECTION_LINE_SEARCH_IMPL_HPP
#define LNOT_BISECTION_LINE_SEARCH_IMPL_HPP

#include <LNOT/LineSearches/BisectionLineSearch.hpp>
#include <LNOT/BasicLinalg/BasicLinalg.hpp>

#include <fmt/core.h>
#include <fmt/format.h>

namespace LNOT
{
	
template<typename T>
BisectionLineSearch<T>::BisectionLineSearch(const Size maxIt, const Scalar firstWolfConditionConst, const Scalar secondWolfConditionConst)
	: Base(maxIt)
	, m_firstWolfConditionConst(firstWolfConditionConst)
	, m_secondWolfConditionConst(secondWolfConditionConst)
{
}

template<typename T>
void BisectionLineSearch<T>::clearWorkSpace()
{
	if (m_newX    != nullptr) { delete[] m_newX;    m_newX    = nullptr; } 
	if (m_newGrad != nullptr) { delete[] m_newGrad; m_newGrad = nullptr; }
	Base::m_workCapacity = 0;
}
	
template<typename T> template<FirstOrderOracle_concept Oracle>
auto BisectionLineSearch<T>::solveImpl(const Scalar* x, const Scalar& fx, const Scalar* gx, const Scalar* s, Oracle& oracle) -> Scalar
{
	constexpr Scalar inf = std::numeric_limits<Scalar>::infinity();
	
	const Size   size     = oracle.getNDims();
	const Scalar sDotGrad = BasicLinalg::inner(s, gx, size);
	
	if (Base::m_workCapacity < size)
	{
		clearWorkSpace();
		Base::m_workCapacity = size;
		m_newX    = new Scalar[Base::m_workCapacity];
		m_newGrad = new Scalar[Base::m_workCapacity];
	}
	
	Scalar alpha = 1;
	Scalar alpha_min = 0;
	Scalar alpha_max = inf;
	
	Base::m_info = Info::FAILURE;
	if (Base::m_out != nullptr) { fmt::print(Base::m_out, "#Bisection LineSearch : \n#Iteration alpha\n"); }
	for (Base::m_nIt=0; Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{		
		if (Base::m_out != nullptr) { fmt::print(Base::m_out, "{} {:10.2e}\n", Base::m_nIt, alpha); }
		if (alpha_max - alpha_min < alpha_max*std::numeric_limits<Scalar>::epsilon()) { break; } // (alpha_min, alpha_max) is empty
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_newX[i] = x[i] + alpha*s[i]; }
		
		oracle.setCurrentPoint(m_newX);
		oracle.getGradient(m_newGrad);
		
		const Scalar fx_new = oracle.getValue();
		const Scalar new_sDotGrad = BasicLinalg::inner(s, m_newGrad, size);
		
		const bool isStepTooShort = new_sDotGrad < m_secondWolfConditionConst*sDotGrad;
		const bool isStepTooLong  = fx_new > fx + m_firstWolfConditionConst*alpha*sDotGrad;
				
		if (not oracle.isFeasible())
		{
			alpha_max = alpha;
			alpha = 0.5*(alpha_max + alpha_min);
		}
		else if (isStepTooLong)
		{
			alpha_max = alpha;
			alpha = 0.5*(alpha_max + alpha_min);
		}
		else if (isStepTooShort)
		{
			alpha_min = alpha;
			alpha = (alpha_max == inf) ? 2*alpha : 0.5*(alpha_max + alpha_min);
		}
		else
		{
			Base::m_info = Info::SUCCESS;
			break;
		}
	}
	return alpha;
}
	
} // namespace LNOT

#endif // LNOT_BISECTION_LINE_SEARCH_IMPL_HPP

