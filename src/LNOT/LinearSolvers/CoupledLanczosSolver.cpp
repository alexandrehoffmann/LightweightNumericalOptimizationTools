#include <LNOT/LinearSolvers/CoupledLanczosSolver.hpp>
#include <LNOT/BasicLinalg/SymmetricDenseMatrixOp.hpp>
#include <LNOT/BasicLinalg/DiagonalPreconditionerOp.hpp>

namespace LNOT
{

template class CoupledLanczosSolver<float>;
template class CoupledLanczosSolver<double>;
template class CoupledLanczosSolver<long double>;

// explicit instanciation for float
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::true_type, Scalar*);
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::true_type, Scalar*);
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::true_type, Scalar*);
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::true_type, Scalar*);

template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::false_type, Scalar*);
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::false_type, Scalar*);
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::false_type, Scalar*);
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::false_type, Scalar*);

template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);

template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);

// explicit instanciation for double
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::true_type, Scalar*);
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::true_type, Scalar*);
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::true_type, Scalar*);
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::true_type, Scalar*);

template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::false_type, Scalar*);
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::false_type, Scalar*);
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::false_type, Scalar*);
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::false_type, Scalar*);

template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);

template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);

// explicit instanciation for long double
template void CoupledLanczosSolver<long double>::solve_impl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::true_type, Scalar*);
template void CoupledLanczosSolver<long double>::solve_impl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::true_type, Scalar*);
template void CoupledLanczosSolver<long double>::solve_impl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::true_type, Scalar*);
template void CoupledLanczosSolver<long double>::solve_impl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::true_type, Scalar*);

template void CoupledLanczosSolver<long double>::solve_impl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::false_type, Scalar*);
template void CoupledLanczosSolver<long double>::solve_impl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::false_type, Scalar*);
template void CoupledLanczosSolver<long double>::solve_impl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::false_type, Scalar*);
template void CoupledLanczosSolver<long double>::solve_impl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::false_type, Scalar*);

template void CoupledLanczosSolver<long double>::solve_impl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void CoupledLanczosSolver<long double>::solve_impl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void CoupledLanczosSolver<long double>::solve_impl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, std::true_type, Scalar*);
template void CoupledLanczosSolver<long double>::solve_impl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, std::true_type, Scalar*);

template void CoupledLanczosSolver<long double>::solve_impl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void CoupledLanczosSolver<long double>::solve_impl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void CoupledLanczosSolver<long double>::solve_impl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, std::false_type, Scalar*);
template void CoupledLanczosSolver<long double>::solve_impl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, std::false_type, Scalar*);

} // namespace LNOT
