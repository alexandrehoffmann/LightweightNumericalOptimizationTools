#ifndef LNOT_DIAGONAL_PRECONDITIONER_OP_HPP
#define LNOT_DIAGONAL_PRECONDITIONER_OP_HPP

#include <concepts>
#include <algorithm>
#include <memory>

namespace LNOT
{

//TODO: add size as a template parameter so it can be a BIC::Fixed
template<std::floating_point Scalar>
class DiagonalPreconditionerOp
{
public:	
	using Size = unsigned int;

	DiagonalPreconditionerOp(const Scalar value, const Size size) : m_invD(std::make_unique<Scalar[]>(size)), m_size(size) { std::fill(m_invD.get(), m_invD.get() + m_size, value); }

	DiagonalPreconditionerOp(const Scalar* A, const Size size) : m_invD(std::make_unique<Scalar[]>(size)), m_size(size) { recompute(A); }
	
	void operator() (const Scalar* x, Scalar* invDx) const;
	
	void recompute(const Scalar* A);
private:
	std::unique_ptr<Scalar[]> m_invD;
	
	Size m_size;
};

} // namespace LNOT

#include <LNOT/BasicLinalg/DiagonalPreconditionerOp_impl.hpp>

#endif // LNOT_DIAGONAL_PRECONDITIONER_OP_HPP
