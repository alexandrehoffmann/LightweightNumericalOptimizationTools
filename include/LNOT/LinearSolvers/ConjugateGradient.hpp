#ifndef LNOT_CONJUGATE_GRADIENT_HPP
#define LNOT_CONJUGATE_GRADIENT_HPP

#include <LNOT/LinearSolvers/LinearSolverBase.hpp>

#include <cmath>

namespace LNOT
{

template<typename T> class ConjugateGradient;

template<typename T> struct LinearSolverTraits< ConjugateGradient<T> > { using Scalar = T;  using Size = unsigned int; };

template<typename T> 
class ConjugateGradient : public LinearSolverBase< ConjugateGradient<T> >
{
	using Base = LinearSolverBase< ConjugateGradient<T> >;
public:
	using Scalar = typename Base::Scalar;
	using Size   = typename Base::Size;
	using Info   = typename Base::Info;

	template<typename HessianOp> using IsHessianOp = typename Base::template IsHessianOp<HessianOp>;
	
	void clearWorkSpace();
	
	template<typename Op, bool solveInPlace> 
	void solve_impl(const Op& H, const Scalar* __restrict__ g, const Size size, std::bool_constant<solveInPlace>, Scalar* __restrict__ x) requires (IsHessianOp<Op>::value);

	Scalar getError        () const { return std::sqrt(m_sqNormR); }
	Scalar getSquaredError () const { return m_sqNormR;            }
private:
	Scalar m_sqNormR = 0;

	Scalar* m_r  = nullptr;
	Scalar* m_p  = nullptr;
	Scalar* m_Hp = nullptr;
};

} // namespace LNOT

#include <LNOT/LinearSolvers/ConjugateGradient_impl.hpp>

#endif // LNOT_CONJUGATE_GRADIENT_HPP
