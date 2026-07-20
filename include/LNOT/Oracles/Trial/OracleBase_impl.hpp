#ifndef LNOT_ORACLE_BASE_IMPL_HPP
#define LNOT_ORACLE_BASE_IMPL_HPP

#include <LNOT/Oracles/OracleBase.hpp>

namespace LNOT
{

template<class Derived>
constexpr auto OracleBase<Derived>::makePrecondOp() -> PrecondOp requires(hasGradient)
{
	if constexpr (hasPrecondOp) 
	{ 
		return CRTP::derived().makePrecondOpImpl(); 
	}
	else
	{
		return PrecondOp(getNDims());
	}
}

template<class Derived>
constexpr auto OracleBase<Derived>::eval(const Scalar* x, Scalar* g) -> Scalar requires(hasGradient)
{
	++m_nGradEval;
	
	if constexpr (detail::CEvalWithGradient<Derived, Scalar>)
	{
		++m_nFuncEval;
		
		return CRTP::derived().evalImpl(x, g); 
	}
	else
	{
		CRTP::derived().getGradient(x, g);
		return eval(x);
	}
}

template<class Derived>
constexpr auto OracleBase<Derived>::eval(const Scalar* x, Scalar* g, PrecondOp& invBop) -> Scalar requires(hasGradient)
{
	if constexpr (hasPrecondOp)
	{
		++m_nPrecondOpEval;
		
		if constexpr (detail::CEvalWithPrecondOp<Derived, Scalar, PrecondOp>)
		{
			++m_nGradEval;
			++m_nFuncEval;
			
			return CRTP::derived().evalImpl(x, g, invBop);
		}
		else
		{
			CRTP::derived().getPrecondOp(x, invBop);
			return eval(x, g);
		}
	}
	else
	{
		return eval(x, g);
	}
}

template<class Derived>
constexpr auto OracleBase<Derived>::eval(const Scalar* x, Scalar* g, HessianOp& Hop) -> Scalar requires(hasHessianOp)
{
	++m_nHessianOpEval;
	
	if constexpr (detail::CEvalWithHessianOp<Derived, Scalar, HessianOp>)
	{
		++m_nGradEval;
		++m_nFuncEval;
		
		return CRTP::derived().evalImpl(x, g, Hop);
	} 
	else
	{
		CRTP::derived().getHessianOp(x, Hop);
		return eval(x, g);
	}
}

template<class Derived>
constexpr auto OracleBase<Derived>::eval(const Scalar* x, Scalar* g, HessianOp& Hop, PrecondOp& invBop) -> Scalar requires(hasHessianOp)
{
	if constexpr (hasPrecondOp)
	{
		if constexpr (detail::CEvalWithHessianAndPrecondOp<Derived, Scalar, HessianOp, PrecondOp>)
		{
			++m_nPrecondOpEval;
			++m_nHessianOpEval;
			++m_nGradEval;
			++m_nFuncEval;
			
			return CRTP::derived().evalImpl(x, g, Hop, invBop);
		}
		else if constexpr (detail::CStandalonePrecondOp<Derived, Scalar, HessianOp>)
		{
			++m_nPrecondOpEval;
			
			CRTP::derived().getPrecondOp(x, invBop);
			return eval(x, g, Hop);
		}
		else
		{
			static_assert(detail::CStandaloneHessianOp<Derived, Scalar, HessianOp>);
			
			++m_nHessianOpEval;
			
			CRTP::derived().getHessianOp(x, invBop);
			return eval(x, g, invBop);
		}
	}
	else
	{
		return eval(x, g, Hop);
	}
}

} // namespace LNOT

#endif // LNOT_ORACLE_BASE_IMPL_HPP
