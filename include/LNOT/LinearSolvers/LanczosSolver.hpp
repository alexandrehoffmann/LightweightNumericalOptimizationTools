#ifndef LNOT_LANCZOS_SOLVER_HPP
#define LNOT_LANCZOS_SOLVER_HPP

#include <LNOT/LinearSolvers/LinearSolverBase.hpp>

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

	Scalar* m_Bz_old = nullptr;
	Scalar* m_Bz     = nullptr;
	Scalar* m_v      = nullptr;
	Scalar* m_p      = nullptr;
	Scalar* m_hat_Bv = nullptr;
};

using LanczosSolverF = LanczosSolver<float>;
using LanczosSolverD = LanczosSolver<double>;

} // namespace LNOT

#include <LNOT/LinearSolvers/LanczosSolver_impl.hpp>

#endif // LNOT_LANCZOS_SOLVER_HPP
