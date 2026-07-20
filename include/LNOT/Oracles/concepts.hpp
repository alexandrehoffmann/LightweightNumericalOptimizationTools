#ifndef LNOT_ORACLE_CONCEPTS_HPP
#define LNOT_ORACLE_CONCEPTS_HPP

#include <concepts>

namespace LNOT
{

template<typename Oracle, typename Scalar, typename Size>
concept COracleInterface = requires(Oracle oracle, const Scalar* x)
{
	{ oracle.getNDimsImpl() }         -> std::convertible_to<Size>;
	{ oracle.setCurrentPointImpl(x) } -> std::same_as<void>;
	{ oracle.isFeasibleImpl() }       -> std::same_as<bool>;
	{ oracle.getValueImpl() }         -> std::convertible_to<Scalar>;
};

template<typename Oracle, typename Scalar, typename Size>
concept CFirstOrderOracleInterface = COracleInterface<Oracle, Scalar, Size> and requires(Oracle oracle, Scalar* g)
{
	{ oracle.getGradientImpl(g) } -> std::same_as<void>;
};

template<typename Oracle, typename Scalar, typename Size>
concept CSecondOrderOracleInterface = CFirstOrderOracleInterface<Oracle, Scalar, Size> and requires(Oracle oracle, const Scalar* d, Scalar* Hd)
{
	{ oracle.getHessianProdImpl(d, Hd) } -> std::same_as<void>;
};

template<typename Oracle, typename Scalar>
concept CPrecondOracleInterface = requires(Oracle oracle, const Scalar* d, Scalar* invBd)
{
	{ oracle.applyPrecondImpl(d, invBd) } -> std::same_as<void>;
};

/////

template<typename Function, typename Scalar>
concept CFunction = std::invocable<Function, const Scalar*> and std::same_as<std::invoke_result_t<Function, const Scalar*>, Scalar>;

template<typename Function, typename Scalar>
concept CGradient = std::invocable<Function, const Scalar*, Scalar*>;

template<typename Function, typename Scalar>
concept CHessianOp = std::invocable<Function, const Scalar*, const Scalar*, Scalar*>;

} // namespace LNOT

#endif // LNOT_ORACLE_CONCEPTS_HPP
