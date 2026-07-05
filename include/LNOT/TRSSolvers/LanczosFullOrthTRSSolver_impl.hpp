#ifndef LNOT_LANCZOS_FULL_ORTH_TRS_SOLVER_IMPL_HPP
#define LNOT_LANCZOS_FULL_ORTH_TRS_SOLVER_IMPL_HPP

#include <LNOT/TRSSolvers/LanczosFullOrthTRSSolver.hpp>
#include <LNOT/BasicLinalg/SymmetricDenseMatrixOp.hpp>
#include <LNOT/BasicLinalg/IdentityPreconditionerOp.hpp>
#include <LNOT/BasicLinalg/DiagonalPreconditionerOp.hpp>

#include <cassert>

namespace LNOT
{

//// explicit template instanciations ////

extern template class LanczosFullOrthTRSSolver<float, 1>;
extern template class LanczosFullOrthTRSSolver<double, 1>;

//// method implementations ////

template<typename T, size_t nPass> 
void LanczosFullOrthTRSSolver<T,nPass>::resizeWorkSpace(const Size newSize)
{
	m_vecSize = newSize;
	m_BV.resize(2*m_vecSize);
	
	if (m_workCapacity < newSize)
	{
		m_workCapacity = newSize;
		
		this->allocateLanczosWorkSpace();
	}
}

template<typename T, size_t nPass> template<typename ASize>
void LanczosFullOrthTRSSolver<T,nPass>::addBvNextImpl(const Scalar scale, const Scalar* unscaledBvNext, const ASize size) requires(isSize<ASize>)
{
	assert(size == m_vecSize);
	
	m_BV.resize(m_BV.size() + m_vecSize);
	
	++m_prev;
	
	Scalar* Bv_curr = getBvCurrImpl();
	
	#pragma omp simd
	for (Size i=0; i!=size; ++i) { Bv_curr[i] = scale*unscaledBvNext[i]; }
}

template<typename T, size_t nPass> template<typename ASize>
void LanczosFullOrthTRSSolver<T,nPass>::reOrthonormalizeImpl(const ASize size, Scalar* unscaledBvNext) requires(isSize<ASize>)
{
	assert(size == m_vecSize);
	
	for (size_t pass=0; pass!=nPass; ++pass)
	{
		for (Size i=0; i!=m_prev+1; ++i)
		{
			const Scalar h = BasicLinalg::inner(getBvImpl(i), unscaledBvNext, size);
			BasicLinalg::axpy(-h, getBvImpl(i), size, unscaledBvNext);
		}
	}
}

} // namespace LNOT

#endif // LNOT_LANCZOS_FULL_ORTH_TRS_SOLVER_IMPL_HPP
