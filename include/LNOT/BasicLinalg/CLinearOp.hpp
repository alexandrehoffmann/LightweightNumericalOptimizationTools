#ifndef LNOT_CLINEAR_OP_HPP
#define LNOT_CLINEAR_OP_HPP

#include <concepts>

namespace LNOT
{
	
template<class Op, class Scalar>
concept CLinearOp = requires(Op Aop, const Scalar* x, Scalar* Ax)
{
	{ Aop.eval(x, Ax) } -> std::same_as<void>;
};
	
} // namespace LNOT

#endif // LNOT_CLINEAR_OP_HPP
