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
	using Self = LanczosTRSSolver<T>;
public:
	LNOT_DEFINE_TRS_SOLVER
	
	LanczosTRSSolver(const Size maxIt = 200000, const Scalar tol = std::numeric_limits<Scalar>::epsilon(), const Size maxItTr = 200000, const Scalar tolTr = std::sqrt(std::numeric_limits<Scalar>::epsilon()));
	
	void clearWorkSpace();
	
	void resizeWorkSpace(const Size newSize);
	
	template<typename HesOp, typename PrecOp, typename ASize> 
	void solveImpl(const HesOp& H, const PrecOp& invB, const Scalar* g, const ASize size, const Scalar& delta, Scalar* x) requires(AreHessianOps<HesOp,PrecOp>::value and IsSize<ASize>::value);
	
	Scalar getErrorImpl        () const { return m_precNormR;             }
	Scalar getSquaredErrorImpl () const { return m_precNormR*m_precNormR; }
	
	Size getMaxItTr() const { return m_maxItTr; }
	
	void setMaxItTr(const Size maxItTr) { m_maxItTr = maxItTr; }
	
	Scalar getLambda() const { return m_lambda; } 
protected:
	LNOT_TRS_SOLVER_ATTRIBUTE
private:
	bool solveBoundary(const Scalar& gamma, const Scalar& delta);

	Size m_maxItTr;

	Scalar m_precNormR  = 0;
	Scalar m_lambda = 0;

	Scalar* m_Bv_old = nullptr;
	Scalar* m_Bv     = nullptr;
	Scalar* m_v      = nullptr;
	Scalar* m_p      = nullptr;
	Scalar* m_Bp     = nullptr;
	Scalar* m_Hp     = nullptr;
	Scalar* m_w      = nullptr;
	
	std::vector<Scalar> m_alpha;
	std::vector<Scalar> m_beta;
	
	std::vector<Scalar> m_invD;
	std::vector<Scalar> m_l;
	std::vector<Scalar> m_h;
};

} // namespace LNOT

#include <LNOT/TRSSolvers/LanczosTRSSolver_impl.hpp>

#endif // LNOT_LANCZOS_TRS_SOLVER_HPP
