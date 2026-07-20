#ifndef LNOT_LINEAR_OP_WRAPPER_HPP
#define LNOT_LINEAR_OP_WRAPPER_HPP

#include <concepts>
#include <algorithm>

namespace LNOT
{

template<class Scalar, class Functor>
class LinearOpWrapper
{
public:
	constexpr LinearOpWrapper(const Functor& functor) : m_op(functor) {}

	constexpr LinearOpWrapper(Functor&& functor) : m_op(std::move(functor)) {}

	constexpr void eval(const Scalar* x, Scalar* precX) { m_op(x, precX); } 
private:	
	Functor m_op;
};

template<class Scalar, std::invocable<const Scalar*, Scalar*> Functor>
LinearOpWrapper<Scalar, std::remove_cvref_t<Functor>> wrappFunctor(Functor&& func) { return LinearOpWrapper<Scalar, std::remove_cvref_t<Functor>>(std::forward<Functor>(func)); }


} // namespace LNOT

#endif // LNOT_LINEAR_OP_WRAPPER_HPP
