#include <LNOT/TRSSolvers/LanczosTRSSolver.hpp>

namespace LNOT
{

template class LanczosTRSSolver<float>;
template class LanczosTRSSolver<double>;
template class LanczosTRSSolver<long double>;

// explicit instanciation for float
template auto LanczosTRSSolver<float>::LanczosTRSSolverBase< LanczosTRSSolver<float> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosTRSSolver<float>::LanczosTRSSolverBase< LanczosTRSSolver<float> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosTRSSolver<float>::LanczosTRSSolverBase< LanczosTRSSolver<float> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosTRSSolver<float>::LanczosTRSSolverBase< LanczosTRSSolver<float> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;

template auto LanczosTRSSolver<float>::LanczosTRSSolverBase< LanczosTRSSolver<float> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosTRSSolver<float>::LanczosTRSSolverBase< LanczosTRSSolver<float> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosTRSSolver<float>::LanczosTRSSolverBase< LanczosTRSSolver<float> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosTRSSolver<float>::LanczosTRSSolverBase< LanczosTRSSolver<float> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;

// explicit instanciation for double
template auto LanczosTRSSolver<double>::LanczosTRSSolverBase< LanczosTRSSolver<double> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosTRSSolver<double>::LanczosTRSSolverBase< LanczosTRSSolver<double> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosTRSSolver<double>::LanczosTRSSolverBase< LanczosTRSSolver<double> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosTRSSolver<double>::LanczosTRSSolverBase< LanczosTRSSolver<double> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;

template auto LanczosTRSSolver<double>::LanczosTRSSolverBase< LanczosTRSSolver<double> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosTRSSolver<double>::LanczosTRSSolverBase< LanczosTRSSolver<double> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosTRSSolver<double>::LanczosTRSSolverBase< LanczosTRSSolver<double> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosTRSSolver<double>::LanczosTRSSolverBase< LanczosTRSSolver<double> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;

} // namespace LNOT
