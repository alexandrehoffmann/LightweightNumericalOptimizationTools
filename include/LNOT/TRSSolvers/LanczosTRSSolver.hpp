#ifndef LNOT_LANCZOS_TRS_SOLVER_HPP
#define LNOT_LANCZOS_TRS_SOLVER_HPP

#include <LNOT/TRSSolvers/TRSSolverBase.hpp>
#include <vector>

namespace LNOT
{

template<typename T> class LanczosTRSSolver;

template<typename T> struct TRSSolverTraits< LanczosTRSSolver<T> > { using Scalar = T; using Size = unsigned int; };

template<typename T>
class LanczosTRSSolver : public TRSSolverBase< LanczosTRSSolver<T> >
{
	using Base = TRSSolverBase< LanczosTRSSolver<T> >;
public:
	using Scalar = typename Base::Scalar;
	using Size   = typename Base::Size;
	using Info   = typename Base::Info;

	template<typename HessianOp> using IsHessianOp = typename Base::template IsHessianOp<HessianOp>;
	
	LanczosTRSSolver(const Size maxIt = 200000, const Scalar tol = std::numeric_limits<Scalar>::epsilon(), const Size maxItTr = 200000, const Scalar tolTr = std::sqrt(std::numeric_limits<Scalar>::epsilon()));
	
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

	Scalar* m_v_old = nullptr;
	Scalar* m_v     = nullptr;
	Scalar* m_p     = nullptr;
	Scalar* m_Hp    = nullptr;
	Scalar* m_w     = nullptr;
	
	std::vector<Scalar> m_alpha;
	std::vector<Scalar> m_beta;
	
	std::vector<Scalar> m_invD;
	std::vector<Scalar> m_l;
	std::vector<Scalar> m_h;
};

} // namespace LNOT

#include <LNOT/TRSSolvers/LanczosTRSSolver_impl.hpp>

#endif // LNOT_LANCZOS_TRS_SOLVER_HPP
