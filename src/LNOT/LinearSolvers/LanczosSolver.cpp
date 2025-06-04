#include <LNOT/LinearSolvers/LanczosSolver.hpp>

namespace LNOT
{

template class LanczosSolver<float>;
template class LanczosSolver<double>;

// explicit instanciation for float
template void LanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);
template void LanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);
template void LanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);
template void LanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);

template void LanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);
template void LanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);
template void LanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);
template void LanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);

// explicit instanciation for double
template void LanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);
template void LanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);
template void LanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);
template void LanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);

template void LanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);
template void LanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);
template void LanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);
template void LanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);

} // namespace LNOT
