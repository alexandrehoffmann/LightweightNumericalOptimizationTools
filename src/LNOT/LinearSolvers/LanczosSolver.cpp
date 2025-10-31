#include <LNOT/LinearSolvers/LanczosSolver.hpp>

namespace LNOT
{

template class LanczosSolver<float>;
template class LanczosSolver<double>;
template class LanczosSolver<long double>;

// explicit instanciation for float
template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);

template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);

template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);

template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);

// explicit instanciation for double
template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);

template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);

template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);

template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);

// explicit instanciation for long double
template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, std::true_type, Scalar*);

template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, std::false_type, Scalar*);

template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, std::true_type, Scalar*);

template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, std::false_type, Scalar*);

} // namespace LNOT
