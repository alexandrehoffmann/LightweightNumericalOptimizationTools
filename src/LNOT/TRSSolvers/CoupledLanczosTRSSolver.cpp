#include <LNOT/TRSSolvers/CoupledLanczosTRSSolver.hpp>

namespace LNOT
{

template class CoupledLanczosTRSSolver<float>;
template class CoupledLanczosTRSSolver<double>;
template class CoupledLanczosTRSSolver<long double>;

// explicit instanciation for float
template void CoupledLanczosTRSSolver<float>::solve(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void CoupledLanczosTRSSolver<float>::solve(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void CoupledLanczosTRSSolver<float>::solve(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void CoupledLanczosTRSSolver<float>::solve(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);

// explicit instanciation for double
template void CoupledLanczosTRSSolver<double>::solve(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void CoupledLanczosTRSSolver<double>::solve(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void CoupledLanczosTRSSolver<double>::solve(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void CoupledLanczosTRSSolver<double>::solve(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
	
} // namespace LNOT
