#ifndef LNOT_SYMMETRIC_DENSE_MATRIX_OP_HPP
#define LNOT_SYMMETRIC_DENSE_MATRIX_OP_HPP

#include <LNOT/BasicLinalg.hpp>

namespace LNOT
{

template<typename Scalar, StorageOrder storageOrder = StorageOrder::ROW_MAJOR, UpLo uplo = UpLo::LOWER>
class SymmetricDenseMatrixOp
{
public:
	using Size = unsigned int;

	SymmetricDenseMatrixOp(const Scalar* A, const Size size) : m_A(A), m_size(size) {}
	
	void operator() (const Scalar* x, Scalar* Ax) const;
	
	constexpr StorageOrder getStorageOrder () const { return storageOrder; }
	constexpr UpLo         getUplo         () const { return uplo;         }
private:
	const Scalar* m_A;
	const Size    m_size;
};

} // namespace LNOT

#include <LNOT/SymmetricDenseMatrixOp_impl.hpp>

#endif // LNOT_SYMMETRIC_DENSE_MATRIX_OP_HPP
