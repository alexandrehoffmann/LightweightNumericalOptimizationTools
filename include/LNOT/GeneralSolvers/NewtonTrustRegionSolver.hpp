#ifndef LNOT_NEWTON_TRUST_REGION_SOLVER_HPP
#define LNOT_NEWTON_TRUST_REGION_SOLVER_HPP

#include <LNOT/GeneralSolvers/SecondOrderSolverBase.hpp>
#include <LNOT/GeneralSolvers/TrustRegionMethodBase.hpp>
#include <LNOT/TRSSolvers/TRSSolverBase.hpp>

namespace LNOT
{
	
template<typename TRSSolver> class NewtonTrustRegionSolver;

template<typename TRSSolver>
struct SecondOrderSolverTraits< NewtonTrustRegionSolver<TRSSolver> >
{
	static_assert(CTRSSolver<TRSSolver>);
	
	using Scalar = typename TRSSolver::Scalar;
	using Size   = typename TRSSolver::Size;
};
	
template<typename TRSSolver>
class NewtonTrustRegionSolver 
	: public SecondOrderSolverBase< NewtonTrustRegionSolver<TRSSolver> >
	, public TrustRegionMethodBase<typename TRSSolver::Scalar, typename TRSSolver::Size>
{
	using Self = NewtonTrustRegionSolver<TRSSolver>;
public:
	LNOT_DEFINE_SECOND_ORDER_SOLVER
	LNOT_DEFINE_TRUST_REGION_SOLVER
	
	void clearWorkSpaceImpl();
	
	template<CSecondOrderOracle Oracle, typename ABool> 
	void solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value);
	
	const TRSSolver& getSubproblemSolver() const { return m_trsSolver; }
	      TRSSolver& getSubproblemSolver()       { return m_trsSolver; }
protected:
	LNOT_SECOND_ORDER_SOLVER_ATTRIBUTE
	LNOT_TRUST_REGION_SOLVER_ATTRIBUTE
private:
	TRSSolver m_trsSolver;
	
	Scalar* m_gk     = nullptr;
	Scalar* m_sk     = nullptr;
	Scalar* m_xTrial = nullptr;
};

} // namespace LNOT

#include <LNOT/GeneralSolvers/NewtonTrustRegionSolver_impl.hpp>

#endif // LNOT_NEWTON_TRUST_REGION_SOLVER_HPP
