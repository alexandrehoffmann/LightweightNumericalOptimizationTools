#ifndef LNOT_DENSE_MATRIX_OP_IMPL_HPP
#define LNOT_DENSE_MATRIX_OP_IMPL_HPP

#include <LNOT/DenseMatrixOp.hpp>

#ifdef LNOT_WITH_BLAS
#include <cblas.h>
#endif // LNOT_WITH_BLAS

namespace LNOT
{

template<typename Scalar, MatrixStorageOrder StorageOrder>
void DenseMatrixOp<Scalar, StorageOrder>::operator() (const Scalar* x, Scalar* Ax) const
{
#ifdef LNOT_WITH_BLAS
	if constexpr (std::is_same<Scalar, float>::value)
	{
		if constexpr (StorageOrder == MatrixStorageOrder::ROW_MAJOR) { cblas_sgemv(CblasRowMajor, CblasNoTrans, blasint(m_nRows), blasint(m_nCols), 1, m_A, blasint(m_nRows), x, 1, 0, Ax, 1); }
		else                                                         { cblas_sgemv(CblasColMajor, CblasNoTrans, blasint(m_nRows), blasint(m_nCols), 1, m_A, blasint(m_nCols), x, 1, 0, Ax, 1); }
	}
	else if constexpr (std::is_same<Scalar, double>::value)
	{
		if constexpr (StorageOrder == MatrixStorageOrder::ROW_MAJOR) { cblas_dgemv(CblasRowMajor, CblasNoTrans, blasint(m_nRows), blasint(m_nCols), 1, m_A, blasint(m_nRows), x, 1, 0, Ax, 1); }
		else                                                         { cblas_dgemv(CblasColMajor, CblasNoTrans, blasint(m_nRows), blasint(m_nCols), 1, m_A, blasint(m_nCols), x, 1, 0, Ax, 1); }
	}
	else
	{
#endif // LNOT_WITH_BLAS
		std::fill(Ax, Ax+m_nRows, 0);
		if constexpr (StorageOrder == MatrixStorageOrder::ROW_MAJOR)
		{
			for (Size i=0; i!=m_nRows; ++i)
			{
				for (Size j=0; j!=m_nCols; ++j)
				{
					Ax[i] += m_A[i*m_nRows + j]*x[j];
				}
			}
		}
		else
		{
			for (Size j=0; j!=m_nCols; ++j)
			{
				for (Size i=0; i!=m_nRows; ++i)
				{
					Ax[i] += m_A[i + j*m_nCols]*x[j];
				}
			}
		}
#ifdef LNOT_WITH_BLAS
	}
#endif // LNOT_WITH_BLAS
}

} // namespace LNOT

#endif // LNOT_DENSE_MATRIX_OP_IMPL_HPP
