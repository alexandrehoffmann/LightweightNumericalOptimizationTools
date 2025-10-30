#ifndef LNOT_IDENTITY_PRECONDITIONER_OP_HPP
#define LNOT_IDENTITY_PRECONDITIONER_OP_HPP

namespace LNOT
{

template<typename Scalar>
class IdentityPreconditionerOp
{
public:
	using Size = unsigned int;
	
	IdentityPreconditionerOp(const Size size) : m_size(size) {}
	
	void operator() (const Scalar* x, Scalar* precX) const { std::copy(x, x + m_size, precX); } 
private:	
	Size m_size;
};

} // namespace LNOT

#endif // LNOT_IDENTITY_PRECONDITIONER_OP_HPP
