#ifndef LNOT_DIAGONAL_PRECONDITIONER_OP_HPP
#define LNOT_DIAGONAL_PRECONDITIONER_OP_HPP

#include <concepts>

namespace LNOT
{

//TODO: add size as a template parameter so it can be a BIC::Fixed
template<std::floating_point Scalar>
class DiagonalPreconditionerOp
{
public:	
	using Size = unsigned int;
	
	DiagonalPreconditionerOp(const Scalar* A, const Size size);
	
	~DiagonalPreconditionerOp() { delete[] m_invD; }
	
	void operator() (const Scalar* x, Scalar* invDx) const;
private:
	Scalar* m_invD;
	
	Size m_size;
};

} // namespace LNOT

#include <LNOT/BasicLinalg/DiagonalPreconditionerOp_impl.hpp>

#endif // LNOT_DIAGONAL_PRECONDITIONER_OP_HPP
