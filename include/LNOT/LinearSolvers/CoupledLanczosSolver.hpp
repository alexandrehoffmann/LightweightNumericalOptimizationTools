#ifndef LNOT_COUPLED_LANCZOS_SOLVER_HPP
#define LNOT_COUPLED_LANCZOS_SOLVER_HPP

#include <LNOT/LinearSolvers/LinearSolverBase.hpp>

namespace LNOT
{

template<typename T> class CoupledLanczosSolver;

template<typename T> struct LinearSolverTraits< CoupledLanczosSolver<T> > { using Scalar = T;  using Size = unsigned int; };

template<typename T> 
class CoupledLanczosSolver : public LinearSolverBase< CoupledLanczosSolver<T> >
{
	using Base = LinearSolverBase< CoupledLanczosSolver<T> >;
public:
	using Scalar = typename Base::Scalar;
	using Size   = typename Base::Size;
	using Info   = typename Base::Info;

	template<typename Hesp>                   using IsHessianOp   = typename Base::template IsHessianOp<Hesp>;
	template<typename HesOp, typename PrecOp> using AreHessianOps = typename Base::template AreHessianOps<HesOp,PrecOp>;
	
	void clearWorkSpace();
	void resizeWorkSpace(const Size newSize);
	
	template<typename Op, bool solveInPlace> 
	void solve_impl(const Op& H, const Scalar* __restrict__ g, const Size size, std::bool_constant<solveInPlace>, Scalar* __restrict__ x) requires (IsHessianOp<Op>::value);
	
	template<typename HesOp, typename PrecOp, bool solveInPlace> 
	void solve_impl(const HesOp& H, const PrecOp& invB, const Scalar* __restrict__ g, const Size size, std::bool_constant<solveInPlace>, Scalar* __restrict__ x) requires (AreHessianOps<HesOp,PrecOp>::value);

	Scalar getError        () const { return m_normR;         }
	Scalar getSquaredError () const { return m_normR*m_normR; }
private:
	Scalar m_normR = 0;

	Scalar* m_v      = nullptr;
	Scalar* m_Bz     = nullptr;
	Scalar* m_q      = nullptr;
	Scalar* m_hat_Bv = nullptr;
};

using CoupledLanczosSolverF = CoupledLanczosSolver<float>;
using CoupledLanczosSolverD = CoupledLanczosSolver<double>;

} // namespace LNOT

#include <LNOT/LinearSolvers/CoupledLanczosSolver_impl.hpp>

#endif // LNOT_COUPLED_LANCZOS_SOLVER_HPP
