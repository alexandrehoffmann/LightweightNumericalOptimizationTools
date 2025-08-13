#ifndef LNOT_SYMMETRIC_DENSE_MATRIX_OP_IMPL_HPP
#define LNOT_SYMMETRIC_DENSE_MATRIX_OP_IMPL_HPP

#include <LNOT/BasicLinalg/SymmetricDenseMatrixOp.hpp>

#ifdef LNOT_WITH_BLAS
#include <cblas.h>
#endif // LNOT_WITH_BLAS

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

} // namespace LNOT

#endif // LNOT_SYMMETRIC_DENSE_MATRIX_OP_IMPL_HPP
