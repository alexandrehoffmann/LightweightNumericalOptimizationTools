#ifndef LNOT_ORACLE_CONCEPTS_HPP
#define LNOT_ORACLE_CONCEPTS_HPP

namespace LNOT
{

template<typename Function, typename Scalar>
concept CFunction = std::invocable<Function, const Scalar*> and std::same_as<std::invoke_result_t<Function, const Scalar*>, Scalar>;

template<typename Function, typename Scalar>
concept CGradient = std::invocable<Function, const Scalar*, Scalar*>;

template<typename Function, typename Scalar>
concept CHessianOp = std::invocable<Function, const Scalar*, const Scalar*, Scalar*>;

} // namespace LNOT

#endif // LNOT_ORACLE_CONCEPTS_HPP
