#ifndef LNOT_SR1_TRUST_REGION_SOLVER_HPP
#define LNOT_SR1_TRUST_REGION_SOLVER_HPP

#include <LNOT/GeneralSolvers/FirstOrderSolverBase.hpp>
#include <LNOT/GeneralSolvers/TrustRegionMethodBase.hpp>
#include <LNOT/TRSSolvers/TRSSolverBase.hpp>

namespace LNOT
{

template<typename TRSSolver> class SR1TrustRegionSolver;

template<typename TRSSolver>
struct FirstOrderSolverTraits< SR1TrustRegionSolver<TRSSolver> >
{
	static_assert(TRSSolver_concept<TRSSolver>);
	
	using Scalar = typename TRSSolver::Scalar;
	using Size   = typename TRSSolver::Size;
};
	
template<typename TRSSolver>
class SR1TrustRegionSolver 
	: public FirstOrderSolverBase< SR1TrustRegionSolver<TRSSolver> >
	, public TrustRegionMethodBase<typename TRSSolver::Scalar, typename TRSSolver::Size>
{
	using Self = SR1TrustRegionSolver<TRSSolver>;
	using TRM  = TrustRegionMethodBase<typename TRSSolver::Scalar, typename TRSSolver::Size>;
public:
	LNOT_DEFINE_FIRST_ORDER_SOLVER
	LNOT_DEFINE_TRUST_REGION_SOLVER
	
	void clearWorkSpaceImpl();
	
	template<FirstOrderOracle_concept Oracle, typename ABool> 
	void solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value);
	
	const TRSSolver& getSubproblemSolver() const { return m_trsSolver; }
	      TRSSolver& getSubproblemSolver()       { return m_trsSolver; }
protected:
	LNOT_FIRST_ORDER_SOLVER_ATTRIBUTE
	LNOT_TRUST_REGION_SOLVER_ATTRIBUTE
private:
	TRSSolver m_trsSolver;
	
	Scalar* m_gk     = nullptr;
	Scalar* m_gkp1   = nullptr;
	Scalar* m_uk     = nullptr;
	Scalar* m_sk     = nullptr;
	Scalar* m_Bk     = nullptr;
	Scalar* m_xTrial = nullptr;
};

} // namespace LNOT

#include <LNOT/GeneralSolvers/SR1TrustRegionSolver_impl.hpp>

#endif // LNOT_SR1_TRUST_REGION_SOLVER_HPP
