#ifndef LNOT_LANCZOS_TRS_SOLVER_IMPL_HPP
#define LNOT_LANCZOS_TRS_SOLVER_IMPL_HPP

#include <LNOT/TRSSolvers/LanczosTRSSolver.hpp>
#include <LNOT/BasicLinalg/SymmetricDenseMatrixOp.hpp>
#include <LNOT/BasicLinalg/IdentityPreconditionerOp.hpp>
#include <LNOT/BasicLinalg/DiagonalPreconditionerOp.hpp>

namespace LNOT
{
	
//// explicit template instanciations ////

extern template class LanczosTRSSolver<float>;
extern template class LanczosTRSSolver<double>;

// explicit instanciation for float
extern template auto LanczosTRSSolver<float>::LanczosTRSSolverBase< LanczosTRSSolver<float> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
extern template auto LanczosTRSSolver<float>::LanczosTRSSolverBase< LanczosTRSSolver<float> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
extern template auto LanczosTRSSolver<float>::LanczosTRSSolverBase< LanczosTRSSolver<float> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
extern template auto LanczosTRSSolver<float>::LanczosTRSSolverBase< LanczosTRSSolver<float> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;

extern template auto LanczosTRSSolver<float>::LanczosTRSSolverBase< LanczosTRSSolver<float> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
extern template auto LanczosTRSSolver<float>::LanczosTRSSolverBase< LanczosTRSSolver<float> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
extern template auto LanczosTRSSolver<float>::LanczosTRSSolverBase< LanczosTRSSolver<float> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
extern template auto LanczosTRSSolver<float>::LanczosTRSSolverBase< LanczosTRSSolver<float> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;

// explicit instanciation for double
extern template auto LanczosTRSSolver<double>::LanczosTRSSolverBase< LanczosTRSSolver<double> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
extern template auto LanczosTRSSolver<double>::LanczosTRSSolverBase< LanczosTRSSolver<double> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
extern template auto LanczosTRSSolver<double>::LanczosTRSSolverBase< LanczosTRSSolver<double> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
extern template auto LanczosTRSSolver<double>::LanczosTRSSolverBase< LanczosTRSSolver<double> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;

extern template auto LanczosTRSSolver<double>::LanczosTRSSolverBase< LanczosTRSSolver<double> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
extern template auto LanczosTRSSolver<double>::LanczosTRSSolverBase< LanczosTRSSolver<double> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
extern template auto LanczosTRSSolver<double>::LanczosTRSSolverBase< LanczosTRSSolver<double> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
extern template auto LanczosTRSSolver<double>::LanczosTRSSolverBase< LanczosTRSSolver<double> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;

//// method implementations ////

template<typename T> 	
void LanczosTRSSolver<T>::clearWorkSpace()
{
	this->clearLanczosWorkSpace();
	if (m_Bv_curr != nullptr) { delete[] m_Bv_curr; m_Bv_curr = nullptr; }
	if (m_Bv_prev != nullptr) { delete[] m_Bv_prev; m_Bv_prev = nullptr; }
}
	
template<typename T> 	
void LanczosTRSSolver<T>::resizeWorkSpace(const Size newSize)
{
	if (m_workCapacity < newSize)
	{
		clearWorkSpace();
		m_workCapacity = newSize;
		
		this->allocateLanczosWorkSpace();
		
		m_Bv_curr = new Scalar[m_workCapacity];
		m_Bv_prev = new Scalar[m_workCapacity];
	}
}

template<typename T> template<typename ASize>
void LanczosTRSSolver<T>::addBvNextImpl(const Scalar scale, const Scalar* unscaledBvNext, const ASize size) requires(IsSize<ASize>::value)
{
	std::move(m_Bv_curr, m_Bv_curr + size, m_Bv_prev);
	
	#pragma omp simd
	for (Size i=0; i!=size; ++i) { m_Bv_curr[i] = scale*unscaledBvNext[i]; }
}
	
} // namespace LNOT

#endif // LNOT_LANCZOS_TRS_SOLVER_IMPL_HPP
