#include <LNOT/TRSSolvers/TruncatedConjugateGradient.hpp>

namespace LNOT
{

template class TruncatedConjugateGradient<float>;
template class TruncatedConjugateGradient<double>;

template void TruncatedConjugateGradient<float>::solve(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void TruncatedConjugateGradient<float>::solve(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void TruncatedConjugateGradient<float>::solve(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void TruncatedConjugateGradient<float>::solve(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);

template void TruncatedConjugateGradient<double>::solve(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void TruncatedConjugateGradient<double>::solve(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void TruncatedConjugateGradient<double>::solve(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void TruncatedConjugateGradient<double>::solve(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
	
} // namespace LNOT
