#ifndef LNOT_LANCZOS_SOLVER_HPP
#define LNOT_LANCZOS_SOLVER_HPP

#include <LNOT/LinearSolvers/LinearSolverBase.hpp>

#include <vector>

namespace LNOT
{

template<typename T> class LanczosSolver;

template<typename T> struct LinearSolverTraits< LanczosSolver<T> > { using Scalar = T;  using Size = unsigned int; };

template<typename T> 
class LanczosSolver : public LinearSolverBase< LanczosSolver<T> >
{
	using Base = LinearSolverBase< LanczosSolver<T> >;
public:
	using Scalar = typename Base::Scalar;
	using Size   = typename Base::Size;
	using Info   = typename Base::Info;

	template<typename HessianOp> using IsHessianOp = typename Base::template IsHessianOp<HessianOp>;
	
	void clearWorkSpace();
	
	template<typename Op, bool solveInPlace> 
	void solve_impl(const Op& H, const Scalar* __restrict__ g, const Size size, std::bool_constant<solveInPlace>, Scalar* __restrict__ x) requires (IsHessianOp<Op>::value);

	Scalar getError        () const { return m_normR;         }
	Scalar getSquaredError () const { return m_normR*m_normR; }
private:
	Scalar m_normR = 0;

	Scalar* m_v_old = nullptr;
	Scalar* m_v     = nullptr;
	Scalar* m_p     = nullptr;
	Scalar* m_w     = nullptr;
};

} // namespace LNOT

#include <LNOT/LinearSolvers/LanczosSolver_impl.hpp>

#endif // LNOT_LANCZOS_SOLVER_HPP
