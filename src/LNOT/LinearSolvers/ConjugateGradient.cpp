#include <LNOT/LinearSolvers/ConjugateGradient.hpp>

namespace LNOT
{

template class ConjugateGradient<float>;
template class ConjugateGradient<double>;

// explicit instanciation for float
template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::true_type, Scalar*);
template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::true_type, Scalar*);
template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::true_type, Scalar*);
template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::true_type, Scalar*);

template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::false_type, Scalar*);
template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::false_type, Scalar*);
template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::false_type, Scalar*);
template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::false_type, Scalar*);

template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);

template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);

// explicit instanciation for double
template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::true_type, Scalar*);
template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::true_type, Scalar*);
template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::true_type, Scalar*);
template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::true_type, Scalar*);

template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::false_type, Scalar*);
template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::false_type, Scalar*);
template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::false_type, Scalar*);
template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::false_type, Scalar*);

template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);

template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);

} // namespace LNOT
