#include <LNOT/TRSSolvers/TruncatedConjugateGradient.hpp>

namespace LNOT
{

//// explicit template instanciations ////

template class TruncatedConjugateGradient<float>;
template class TruncatedConjugateGradient<double>;

// explicit instanciation for float
template void TruncatedConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void TruncatedConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void TruncatedConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void TruncatedConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);

template void TruncatedConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void TruncatedConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void TruncatedConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void TruncatedConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);

// explicit instanciation for double
template void TruncatedConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void TruncatedConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void TruncatedConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void TruncatedConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);

template void TruncatedConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void TruncatedConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void TruncatedConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
template void TruncatedConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);

} // namespace LNOT
