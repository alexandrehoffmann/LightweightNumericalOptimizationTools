#ifndef LNOT_DENSE_MATRIX_OP_HPP
#define LNOT_DENSE_MATRIX_OP_HPP

#include <algorithm>

namespace LNOT
{

enum class MatrixStorageOrder { ROW_MAJOR, COL_MAJOR };

template<typename Scalar, MatrixStorageOrder StorageOrder = MatrixStorageOrder::ROW_MAJOR>
class DenseMatrixOp
{
public:
	using Size = unsigned int;

	DenseMatrixOp(const Scalar* A, const Size nRows, const Size nCols) : m_A(A), m_nRows(nRows), m_nCols(nCols) {}
	
	void operator() (const Scalar* x, Scalar* Ax) const;
private:
	const Scalar* m_A;
	const Size    m_nRows;
	const Size    m_nCols;
};

} // namespace LNOT

#include <LNOT/DenseMatrixOp_impl.hpp>

#endif // LNOT_DENSE_MATRIX_OP_HPP
