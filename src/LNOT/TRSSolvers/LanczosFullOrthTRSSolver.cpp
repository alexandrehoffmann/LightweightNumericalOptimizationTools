#include <LNOT/TRSSolvers/LanczosFullOrthTRSSolver.hpp>

namespace LNOT
{


template class LanczosFullOrthTRSSolver<float, 1>;
template class LanczosFullOrthTRSSolver<double, 1>;

// explicit instanciation for float
template auto LanczosFullOrthTRSSolver<float, 1>::LanczosTRSSolverBase< LanczosFullOrthTRSSolver<float, 1> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosFullOrthTRSSolver<float, 1>::LanczosTRSSolverBase< LanczosFullOrthTRSSolver<float, 1> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosFullOrthTRSSolver<float, 1>::LanczosTRSSolverBase< LanczosFullOrthTRSSolver<float, 1> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosFullOrthTRSSolver<float, 1>::LanczosTRSSolverBase< LanczosFullOrthTRSSolver<float, 1> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;

template auto LanczosFullOrthTRSSolver<float, 1>::LanczosTRSSolverBase< LanczosFullOrthTRSSolver<float, 1> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosFullOrthTRSSolver<float, 1>::LanczosTRSSolverBase< LanczosFullOrthTRSSolver<float, 1> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosFullOrthTRSSolver<float, 1>::LanczosTRSSolverBase< LanczosFullOrthTRSSolver<float, 1> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosFullOrthTRSSolver<float, 1>::LanczosTRSSolverBase< LanczosFullOrthTRSSolver<float, 1> >::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;

// explicit instanciation for double
template auto LanczosFullOrthTRSSolver<double, 1>::LanczosTRSSolverBase< LanczosFullOrthTRSSolver<double, 1> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosFullOrthTRSSolver<double, 1>::LanczosTRSSolverBase< LanczosFullOrthTRSSolver<double, 1> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosFullOrthTRSSolver<double, 1>::LanczosTRSSolverBase< LanczosFullOrthTRSSolver<double, 1> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosFullOrthTRSSolver<double, 1>::LanczosTRSSolverBase< LanczosFullOrthTRSSolver<double, 1> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;

template auto LanczosFullOrthTRSSolver<double, 1>::LanczosTRSSolverBase< LanczosFullOrthTRSSolver<double, 1> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosFullOrthTRSSolver<double, 1>::LanczosTRSSolverBase< LanczosFullOrthTRSSolver<double, 1> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosFullOrthTRSSolver<double, 1>::LanczosTRSSolverBase< LanczosFullOrthTRSSolver<double, 1> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;
template auto LanczosFullOrthTRSSolver<double, 1>::LanczosTRSSolverBase< LanczosFullOrthTRSSolver<double, 1> >::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) -> Scalar;

} // namespace LNOT
