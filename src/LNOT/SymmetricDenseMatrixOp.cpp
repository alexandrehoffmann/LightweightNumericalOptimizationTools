#include <LNOT/SymmetricDenseMatrixOp.hpp>

namespace LNOT
{

template class SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>;
template class SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>;
template class SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>;
template class SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>;

template class SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>;
template class SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>;
template class SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>;
template class SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>;

} // namespace LNOT
