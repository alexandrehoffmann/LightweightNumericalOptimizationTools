#ifndef LNOT_BACKTRACKING_LINESEARCH_IMPL_HPP
#define LNOT_BACKTRACKING_LINESEARCH_IMPL_HPP

#include <LNOT/LineSearches/BacktrackingLineSearch.hpp>
#include <LNOT/BasicLinalg/BasicLinalg.hpp>

#include <fmt/core.h>
#include <fmt/format.h>

namespace LNOT
{

template<typename T>
void BacktrackingLineSearch<T>::clearWorkSpace()
{
	if (m_newX != nullptr) { delete[] m_newX; m_newX = nullptr; }
	Base::m_workCapacity = 0;
}

template<typename T> template<FirstOrderOracle_concept Oracle>
auto BacktrackingLineSearch<T>::solve(const Scalar* x, const Scalar fx, const Scalar* gx, const Scalar* s, Oracle& oracle) -> Scalar
{	
	const Size   size     = oracle.getNDims();
	const Scalar sDotGrad = BasicLinalg::inner(s, gx, size);
	
	if (Base::m_workCapacity < size)
	{
		clearWorkSpace();
		Base::m_workCapacity = size;
		m_newX = new Scalar[Base::m_workCapacity];
	}
	
	Scalar alpha = 1;
	
	Base::m_info = Info::FAILURE;
	if (Base::m_out != nullptr) { fmt::print(Base::m_out, "#Backtracking LineSearch : \n#Iteration alpha f(x+alpha s) f(x) tol\n"); }
	for (Base::m_nIt=0; Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{
		if (alpha < std::numeric_limits<Scalar>::epsilon()) { break; }
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_newX[i] = x[i] + alpha*s[i]; }
		
		oracle.setCurrentPoint(m_newX);
		const Scalar fx_new = oracle.getValue();
		
		if (Base::m_out != nullptr) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e} {:10.2e} {:10.2e}\n", Base::m_nIt, alpha, fx_new, fx, m_armijoConditionConst*alpha*sDotGrad); }
		
		if (not oracle.isFeasible())
		{
			alpha *= m_tau;
		}
		else if (fx_new - fx > m_armijoConditionConst*alpha*sDotGrad)
		{
			alpha *= m_tau;
		}
		else
		{
			Base::m_info = Info::SUCCESS;
			break;
		}
	}
	return alpha;
}
	
	
}

#endif // LNOT_BACKTRACKING_LINESEARCH_IMPL_HPP

