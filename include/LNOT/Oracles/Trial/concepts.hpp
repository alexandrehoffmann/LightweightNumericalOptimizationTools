#ifndef LNOT_ORACLE_CONCEPTS_HPP
#define LNOT_ORACLE_CONCEPTS_HPP

#include <concepts>

namespace LNOT
{
	
namespace detail
{
	template<typename Oracle, typename Scalar>
	concept CEvalWithGradient = requires(Oracle oracle, const Scalar* x, Scalar* g) 
	{ 
		{ oracle.evalImpl(x, g) } -> std::same_as<Scalar>; 
	};
	
	template<typename Oracle, typename Scalar>
	concept CStandaloneGradient = requires(Oracle oracle, const Scalar* x, Scalar* g) 
	{ 
		{ oracle.getGradient(x, g) } -> std::same_as<void>;
	};
	
	template<typename Oracle, typename Scalar, typename HessianOp>
	concept CEvalWithHessianOp = requires(Oracle oracle, const Scalar* x, Scalar* g, HessianOp& Hop) 
	{ 
		{ oracle.evalImpl(x, g, Hop) } -> std::same_as<Scalar>; 
	};
	
	template<typename Oracle, typename Scalar, typename HessianOp>
	concept CStandaloneHessianOp = requires(Oracle oracle, const Scalar* x, HessianOp& Hop) 
	{ 
		oracle.getHessianOp(x, Hop); 
	};
	
	template<typename Oracle, typename Scalar, typename HessianOp, typename PrecondOp>
	concept CEvalWithHessianAndPrecondOp = requires(Oracle oracle, const Scalar* x, Scalar* g, HessianOp& Hop, PrecondOp& invBop) 
	{ 
		{ oracle.evalImpl(x, g, Hop, invBop) } -> std::same_as<Scalar>;
	};
	
	template<typename Oracle, typename Scalar, typename PrecondOp>
	concept CEvalWithPrecondOp = requires(Oracle oracle, const Scalar* x, Scalar* g, PrecondOp& invBop) 
	{ 
		{ oracle.evalImpl(x, g, invBop) } -> std::same_as<Scalar>; 
	};
	
	template<typename Oracle, typename Scalar, typename PrecondOp>
	concept CStandalonePrecondOp = requires(Oracle oracle, const Scalar* x, PrecondOp& invBop) 
	{ 
		{ oracle.getPrecondOp(x, invBop) } -> std::same_as<void>; 
	};
	
} // namespace detail

template<typename Oracle, typename Scalar>
concept CFirstOrderInterface = detail::CEvalWithGradient<Oracle, Scalar> or detail::CStandaloneGradient<Oracle, Scalar>;

template<typename Oracle, typename Scalar, typename PrecondOp>
concept CPrecondFirstOrderInterface = (CFirstOrderInterface<Oracle, Scalar> and detail::CStandalonePrecondOp<Oracle, Scalar, PrecondOp>) or detail::CEvalWithPrecondOp<Oracle, Scalar, PrecondOp>;

template<typename Oracle, typename Scalar, typename HessianOp>
concept CSecondOrderInterface = (CFirstOrderInterface<Oracle, Scalar> and detail::CStandaloneHessianOp<Oracle, Scalar, HessianOp>) or detail::CEvalWithHessianOp<Oracle, Scalar, HessianOp>;

template<typename Oracle, typename Scalar, typename HessianOp, typename PrecondOp>
concept CPrecondhSecondOrderInterface = (CSecondOrderInterface<Oracle, Scalar, HessianOp> and detail::CStandalonePrecondOp<Oracle, Scalar, PrecondOp>) or (CPrecondFirstOrderInterface<Oracle, Scalar, PrecondOp> and detail::CStandaloneHessianOp<Oracle, Scalar, HessianOp>) or detail::CEvalWithHessianAndPrecondOp<Oracle, Scalar, HessianOp, PrecondOp>; 

/////

template<typename Function, typename Scalar>
concept CFunction = std::invocable<Function, const Scalar*> and std::same_as<std::invoke_result_t<Function, const Scalar*>, Scalar>;

template<typename Function, typename Scalar>
concept CGradient = std::invocable<Function, const Scalar*, Scalar*>;

template<typename Function, typename Scalar>
concept CHessianOp = std::invocable<Function, const Scalar*, const Scalar*, Scalar*>;

} // namespace LNOT

#endif // LNOT_ORACLE_CONCEPTS_HPP
