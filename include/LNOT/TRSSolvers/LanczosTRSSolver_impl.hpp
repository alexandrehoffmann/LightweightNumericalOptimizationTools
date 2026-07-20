#ifndef LNOT_LANCZOS_TRS_SOLVER_IMPL_HPP
#define LNOT_LANCZOS_TRS_SOLVER_IMPL_HPP

#include <LNOT/TRSSolvers/LanczosTRSSolver.hpp>
namespace LNOT
{
	
//// explicit template instanciations ////

extern template class LanczosTRSSolver<float>;
extern template class LanczosTRSSolver<double>;

//// method implementations ////

template<typename T> 	
void LanczosTRSSolver<T>::resizeWorkSpace(const Size newSize)
{
	if (not m_Bv_curr or m_workCapacity < newSize)
	{
		m_workCapacity = newSize;
		
		this->allocateLanczosWorkSpace();
		
		m_Bv_curr = std::make_unique<Scalar[]>(m_workCapacity);
		m_Bv_prev = std::make_unique<Scalar[]>(m_workCapacity);
	}
}

template<typename T> template<typename ASize>
void LanczosTRSSolver<T>::addBvNextImpl(const Scalar scale, const Scalar* unscaledBvNext, const ASize size) requires(isSize<ASize>)
{
	std::move(m_Bv_curr.get(), m_Bv_curr.get() + size, m_Bv_prev.get());
	
	#pragma omp simd
	for (Size i=0; i!=size; ++i) { m_Bv_curr[i] = scale*unscaledBvNext[i]; }
}
	
} // namespace LNOT

#endif // LNOT_LANCZOS_TRS_SOLVER_IMPL_HPP
