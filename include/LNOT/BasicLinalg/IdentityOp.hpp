#ifndef LNOT_IDENTITY_OP_HPP
#define LNOT_IDENTITY_OP_HPP

#include <concepts>
#include <algorithm>

namespace LNOT
{

class IdentityOp
{
public:
	using Size = unsigned int;
	
	constexpr IdentityOp(const Size size) : m_size(size) {}
	
	template<class Scalar>
	constexpr void eval(const Scalar* x, Scalar* precX) const { std::copy(x, x + m_size, precX); } 
private:	
	Size m_size;
};

} // namespace LNOT

#endif // LNOT_IDENTITY_OP_HPP
