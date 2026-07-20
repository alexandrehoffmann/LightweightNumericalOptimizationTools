#ifndef LNOT_SYMMETRIC_DENSE_MATRIX_OP_IMPL_HPP
#define LNOT_SYMMETRIC_DENSE_MATRIX_OP_IMPL_HPP

#include <LNOT/BasicLinalg/SymmetricDenseMatrixOp.hpp>

namespace LNOT
{
	
//// explicit template instanciations ////

extern template class SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>; 
extern template class SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>; 
extern template class SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>; 
extern template class SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>; 

extern template class SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>; 
extern template class SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>; 
extern template class SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>; 
extern template class SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>; 

extern template class SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::LOWER>; 
extern template class SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::UPPER>; 
extern template class SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::LOWER>; 
extern template class SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::UPPER>; 

//// method implementation ////

template<class Scalar, StorageOrder storageOrder, UpLo uplo>
constexpr auto SymmetricDenseMatrixOp<Scalar, storageOrder, uplo>::getFlatIndex(Size i, Size j) const -> Size
{
	if ((uplo == UpLo::UPPER and i < j) or (uplo == UpLo::LOWER and i > j))
	{
		std::swap(i, j);
	}
	
	const Size iStride = storageOrder == StorageOrder::ROW_MAJOR ? m_size  : Size(1);
	const Size jStride = storageOrder == StorageOrder::ROW_MAJOR ? Size(1) : m_size;
	
	return iStride*i + jStride*j;
}

} // namespace LNOT

#endif // LNOT_SYMMETRIC_DENSE_MATRIX_OP_IMPL_HPP
