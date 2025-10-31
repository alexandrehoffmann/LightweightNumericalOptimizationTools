#ifndef LNOT_ORACLE_WRAPPER_IMPL_HPP
#define LNOT_ORACLE_WRAPPER_IMPL_HPP

#include <LNOT/Oracles/OracleWrapper.hpp>

namespace LNOT
{

template<typename T, typename S, class Function, class Gradient, class HessianOp, class PrecondOp> 
void OracleWrapper<T, S, Function, Gradient, HessianOp, PrecondOp>::applyPrecondImpl(const Scalar* d, Scalar* invBd) const requires (hasHessianProd) 
{ 
	if constexpr (hasApplyPrecond) { m_precondOp(m_x, d, invBd);           } 
	else                           { Base::applyIdentityPrecond(d, invBd); } 
}
	
} // namespace LNOT

#endif // LNOT_ORACLE_WRAPPER_IMPL_HPP
