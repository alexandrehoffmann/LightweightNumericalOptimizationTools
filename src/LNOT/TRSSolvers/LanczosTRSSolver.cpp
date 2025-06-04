#include <LNOT/TRSSolvers/LanczosTRSSolver.hpp>

namespace LNOT
{

template class LanczosTRSSolver<float>;
template class LanczosTRSSolver<double>;

template void LanczosTRSSolver<float>::solve(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void LanczosTRSSolver<float>::solve(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void LanczosTRSSolver<float>::solve(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void LanczosTRSSolver<float>::solve(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);

template void LanczosTRSSolver<double>::solve(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void LanczosTRSSolver<double>::solve(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void LanczosTRSSolver<double>::solve(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void LanczosTRSSolver<double>::solve(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
	
} // namespace LNOT
