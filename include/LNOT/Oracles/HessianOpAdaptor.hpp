#ifndef LNOT_HESSIAN_OP_ADAPTOR_HPP
#define LNOT_HESSIAN_OP_ADAPTOR_HPP

#include <LNOT/Oracles/concepts.hpp>

namespace LNOT
{
	
template<class Derived> struct OracleTraits;

template<class Oracle>
class HessianOpAdaptor
{
public:
	using Scalar = typename OracleTraits<Oracle>::Scalar;

	constexpr HessianOpAdaptor(Oracle& oracle) : m_oracle(oracle) { m_oracle.initHessianProd(); }
	
	constexpr void eval(const Scalar* d, Scalar* Hd) { m_oracle.getHessianProdImpl(d, Hd); }
private:
	Oracle& m_oracle;
}; 

} // namespace LNOT

#endif // LNOT_HESSIAN_OP_ADAPTOR_HPP
