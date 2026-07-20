#ifndef LNOT_SYMMETRIC_DENSE_MATRIX_OP_HPP
#define LNOT_SYMMETRIC_DENSE_MATRIX_OP_HPP

#include <LNOT/BasicLinalg/BasicLinalg.hpp>
#include <BIC/Core.hpp>

#include <type_traits>
#include <span>

namespace LNOT
{
	
template<class T, StorageOrder storageOrder = StorageOrder::ROW_MAJOR, UpLo uplo = UpLo::LOWER>
class SymmetricDenseMatrixOp
{
public:
    using Element         = T;
    using const_Element   = std::add_const_t<Element>;
    using Scalar          = std::remove_cv_t<Element>;
    using Pointer         = Element*;
    using Reference       = Element&;
    using const_Reference = std::add_const_t<Element>&;

    static constexpr bool isConst = std::same_as<Element, const_Element>;
	
	using Size = unsigned int;

	constexpr SymmetricDenseMatrixOp(Pointer A, const Size size) : m_A(A), m_size(size) {}
	
	constexpr void eval(const Scalar* x, Scalar* Ax) const { BasicLinalg::symMatrixVectorProd(storageOrder, uplo, Scalar(1), m_A, x, m_size, BIC::fixed<bool,false>, Ax); }
	
	constexpr Size getSize() const { return m_size; }

	constexpr       Reference operator()(const Size i, const Size j) requires(not isConst) { return m_A[getFlatIndex(i, j)]; }
	constexpr const_Reference operator()(const Size i, const Size j) const                 { return m_A[getFlatIndex(i, j)]; }

	constexpr std::span<      Element> getElements() requires(not isConst) { return std::span(m_A, m_size*m_size); }
	constexpr std::span<const_Element> getElements() const                 { return std::span(m_A, m_size*m_size); }
	
	constexpr StorageOrder getStorageOrder () const { return storageOrder; }
	constexpr UpLo         getUplo         () const { return uplo;         }
private:
	constexpr Size getFlatIndex(Size i, Size j) const;

	Pointer m_A;
	
	Size m_size;
};

} // namespace LNOT

#include <LNOT/BasicLinalg/SymmetricDenseMatrixOp_impl.hpp>

#endif // LNOT_SYMMETRIC_DENSE_MATRIX_OP_HPP
