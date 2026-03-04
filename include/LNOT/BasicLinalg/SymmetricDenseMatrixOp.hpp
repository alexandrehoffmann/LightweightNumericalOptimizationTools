#ifndef LNOT_SYMMETRIC_DENSE_MATRIX_OP_HPP
#define LNOT_SYMMETRIC_DENSE_MATRIX_OP_HPP

#include <LNOT/BasicLinalg/BasicLinalg.hpp>
#include <BIC/Core.hpp>

#include <concepts>

namespace LNOT
{

//TODO: add size as a template parameter so it can be a BIC::Fixed
template<std::floating_point Scalar, StorageOrder storageOrder = StorageOrder::ROW_MAJOR, UpLo uplo = UpLo::LOWER>
class SymmetricDenseMatrixOp
{
public:
	using Size = unsigned int;

	SymmetricDenseMatrixOp(const Scalar* A, const Size size) : m_A(A), m_size(size) {}
	
	void operator() (const Scalar* x, Scalar* Ax) const { BasicLinalg::symMatrixVectorProd(storageOrder, uplo, Scalar(1), m_A, x, m_size, BIC::fixed<bool,false>, Ax); }
	
	constexpr StorageOrder getStorageOrder () const { return storageOrder; }
	constexpr UpLo         getUplo         () const { return uplo;         }
private:
	const Scalar* m_A;
	
	Size m_size;
};

} // namespace LNOT

#include <LNOT/BasicLinalg/SymmetricDenseMatrixOp_impl.hpp>

#endif // LNOT_SYMMETRIC_DENSE_MATRIX_OP_HPP
