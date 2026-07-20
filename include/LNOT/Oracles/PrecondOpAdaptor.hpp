#ifndef LNOT_PRECOND_OP_ADAPTOR_HPP
#define LNOT_PRECOND_OP_ADAPTOR_HPP

#include <LNOT/Oracles/concepts.hpp>

namespace LNOT
{
	
template<class Derived> struct OracleTraits;

template<class Oracle>
class PrecondOpAdaptor
{
public:
	using Scalar = typename OracleTraits<Oracle>::Scalar;

	constexpr PrecondOpAdaptor(Oracle& oracle) : m_oracle(oracle) { m_oracle.initApplyPrecond(); }
	
	constexpr void eval(const Scalar* d, Scalar* invBd) { m_oracle.applyPrecondImpl(d, invBd); }
private:
	Oracle& m_oracle;
}; 

} // namespace LNOT

#endif // LNOT_PRECOND_OP_ADAPTOR_HPP
