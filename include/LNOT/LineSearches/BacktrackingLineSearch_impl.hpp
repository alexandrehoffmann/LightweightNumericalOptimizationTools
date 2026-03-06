#ifndef LNOT_BACKTRACKING_LINESEARCH_IMPL_HPP
#define LNOT_BACKTRACKING_LINESEARCH_IMPL_HPP

#include <LNOT/LineSearches/BacktrackingLineSearch.hpp>
#include <LNOT/BasicLinalg/BasicLinalg.hpp>
#include <LNOT/FloatingPoint/FPComparator.hpp>

#include <fmt/core.h>
#include <fmt/format.h>

namespace LNOT
{

//// explicit template instanciations ////

extern template class BacktrackingLineSearch<float>;
extern template class BacktrackingLineSearch<double>;
extern template class BacktrackingLineSearch<long double>;

//// method implementations ////

template<typename T>
void BacktrackingLineSearch<T>::clearWorkSpace()
{
	if (m_newX != nullptr) { delete[] m_newX; m_newX = nullptr; }
	m_workCapacity = 0;
}

template<typename T> template<CFirstOrderOracle Oracle>
auto BacktrackingLineSearch<T>::solveImpl(const Scalar* x, const Scalar fx, const Scalar* gx, const Scalar* s, Oracle& oracle) -> Scalar
{	
	using Oracle_Size = typename Oracle::Size;
	
	const Oracle_Size size     = oracle.getNDims();
	const Scalar      sDotGrad = BasicLinalg::inner(s, gx, size);
	
	const FPComparator<Scalar> cmp;
	
	if (m_workCapacity < size)
	{
		clearWorkSpace();
		m_workCapacity = size;
		m_newX = new Scalar[m_workCapacity];
	}
	
	const Scalar normS = BasicLinalg::norm(s, size);
	
	Scalar alpha(1);
	
	m_info = Info::FAILURE;
	if (m_out != nullptr) { fmt::println(m_out, "#Backtracking LineSearch : \n#Iteration alpha f(x+alpha s) f(x) tol"); }
	for (m_nIt=0; m_nIt!=m_maxIt; ++m_nIt)
	{
		if (not cmp.isDefPositive(alpha*normS)) { break; }
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_newX[i] = x[i] + alpha*s[i]; }
		
		oracle.setCurrentPoint(m_newX);
		const Scalar fx_new = oracle.getValue();
		
		if (m_out != nullptr) { fmt::println(m_out, "{} {:10.2e} {:10.2e} {:10.2e} {:10.2e}", m_nIt, alpha, fx_new, fx, m_armijoConditionConst*alpha*sDotGrad); std::fflush(m_out); }
		
		if (not oracle.isFeasible())
		{
			alpha *= m_tau;
		}
		else if (cmp.isDefGreaterThan(fx_new - fx, m_armijoConditionConst*alpha*sDotGrad))
		{
			alpha *= m_tau;
		}
		else
		{
			m_info = Info::SUCCESS;
			break;
		}
	}
	return alpha;
}
	
	
}

#endif // LNOT_BACKTRACKING_LINESEARCH_IMPL_HPP

