#ifndef LNOT_COUPLED_LANCZOS_TRS_SOLVER_HPP
#define LNOT_COUPLED_LANCZOS_TRS_SOLVER_HPP

#include <LNOT/TRSSolvers/TRSSolverBase.hpp>
#include <vector>

namespace LNOT
{

template<typename T> class CoupledLanczosTRSSolver;

template<typename T> struct TRSSolverTraits< CoupledLanczosTRSSolver<T> > { using Scalar = T; using Size = unsigned int; };

template<typename T>
class CoupledLanczosTRSSolver : public TRSSolverBase< CoupledLanczosTRSSolver<T> >
{
	using Base = TRSSolverBase< CoupledLanczosTRSSolver<T> >;
public:
	using Scalar = typename Base::Scalar;
	using Size   = typename Base::Size;
	using Info   = typename Base::Info;

	template<typename HessianOp> using IsHessianOp = typename Base::template IsHessianOp<HessianOp>;
	
	CoupledLanczosTRSSolver(const Size maxIt = 200000, const Scalar tol = std::numeric_limits<Scalar>::epsilon(), const Size maxItTr = 200000, const Scalar tolTr = std::sqrt(std::numeric_limits<Scalar>::epsilon()));
	
	void clearWorkSpace();
	
	template<typename Op> 
	void solve(const Op& H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) requires (IsHessianOp<Op>::value);
	
	Scalar getError        () const { return m_normR;         }
	Scalar getSquaredError () const { return m_normR*m_normR; }
	
	Size getMaxItTr() const { return m_maxItTr; }
	
	void setMaxItTr(const Size maxItTr) { m_maxItTr = maxItTr; }
	
	Scalar getLambda() const { return m_lambda; } 
private:
	bool solveBoundary(const Scalar& __restrict__ gamma, const Scalar& __restrict__ delta);

	Size m_maxItTr;

	Scalar m_normR  = 0;
	Scalar m_lambda = 0;
	
	Scalar* m_v     = nullptr;
	Scalar* m_q     = nullptr;
	Scalar* m_w     = nullptr;
	Scalar* m_w_old = nullptr; // only used to compute model reduction
	
	std::vector<Scalar> m_alpha;
	std::vector<Scalar> m_beta;
	
	std::vector<Scalar> m_invD;
	std::vector<Scalar> m_l;
	std::vector<Scalar> m_h;
};

} // namespace LNOT

#include <LNOT/TRSSolvers/CoupledLanczosTRSSolver_impl.hpp>

#endif // LNOT_COUPLED_LANCZOS_TRS_SOLVER_HPP
