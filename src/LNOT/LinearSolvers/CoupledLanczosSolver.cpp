#include <LNOT/LinearSolvers/CoupledLanczosSolver.hpp>
#include <LNOT/SymmetricDenseMatrixOp.hpp>

namespace LNOT
{

template class CoupledLanczosSolver<float>;
template class CoupledLanczosSolver<double>;

// explicit instanciation for float
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);

template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);
template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);

// explicit instanciation for double
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);

template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);
template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);

} // namespace LNOT
