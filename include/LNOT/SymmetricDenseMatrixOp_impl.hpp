#ifndef LNOT_SYMMETRIC_DENSE_MATRIX_OP_IMPL_HPP
#define LNOT_SYMMETRIC_DENSE_MATRIX_OP_IMPL_HPP

#include <LNOT/SymmetricDenseMatrixOp.hpp>

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

//// method implementations ////

template<typename Scalar, StorageOrder storageOrder, UpLo uplo>
void SymmetricDenseMatrixOp<Scalar, storageOrder, uplo>::operator() (const Scalar* __restrict__ x, Scalar* __restrict__ Ax) const
{
#ifdef LNOT_WITH_BLAS
	if constexpr      (std::is_same<Scalar, float>::value)  { cblas_ssymv(CBLAS_ORDER(storageOrder), CBLAS_UPLO(uplo), blasint(m_size), 1, m_A, blasint(m_size), x, 1, 0, Ax, 1); }
	else if constexpr (std::is_same<Scalar, double>::value) { cblas_dsymv(CBLAS_ORDER(storageOrder), CBLAS_UPLO(uplo), blasint(m_size), 1, m_A, blasint(m_size), x, 1, 0, Ax, 1); }
	else
	{
#endif // LNOT_WITH_BLAS
		std::fill(Ax, Ax+m_size, 0);
		
		// if uplo == UPPER -> we need to transpose, i.e. switch the strides
		const Size iStride = storageOrder == StorageOrder::ROW_MAJOR and uplo == UpLo::LOWER ? m_size : 1;
		const Size jStride = storageOrder == StorageOrder::ROW_MAJOR and uplo == UpLo::LOWER ?      1 : m_size;
		
		for (Size i=0; i!=m_size; ++i)
		{
			for (Size j=0; j!=i; ++j)
			{
				Ax[i] += m_A[i*iStride + j*jStride]*x[j];
			}
			Ax[i] += m_A[i*iStride + i*jStride]*x[i];
			for (Size j=i+1; j!=m_size; ++j)
			{
				Ax[i] += m_A[j*iStride + i*jStride]*x[j];
			}
		}
#ifdef LNOT_WITH_BLAS
	}
#endif // LNOT_WITH_BLAS
}

} // namespace LNOT

#endif // LNOT_SYMMETRIC_DENSE_MATRIX_OP_IMPL_HPP
