#ifndef LNOT_LSR1_TRUST_REGION_SOLVER_HPP
#define LNOT_LSR1_TRUST_REGION_SOLVER_HPP

#include <LNOT/GeneralSolvers/FirstOrderSolverBase.hpp>
#include <LNOT/GeneralSolvers/TrustRegionMethodBase.hpp>
#include <LNOT/GeneralSolvers/LimitedMemorySolverBase.hpp>
#include <LNOT/TRSSolvers/TRSSolverBase.hpp>
#include <LNOT/ConvergenceCriterions/L2Norm.hpp>

#include <memory>

namespace LNOT
{

template<typename TRSSolver, typename ConvergenceCriterion> class LSR1TrustRegionSolver;

template<typename TRSSolver, typename ConvergenceCriterion>
struct FirstOrderSolverTraits< LSR1TrustRegionSolver<TRSSolver, ConvergenceCriterion> >
{
	static_assert(CTRSSolver<TRSSolver>);
	static_assert(CConvergenceCriterion<ConvergenceCriterion>);
	
	static_assert(std::is_same<typename TRSSolver::Scalar, typename ConvergenceCriterion::Scalar>::value);
	static_assert(std::is_same<typename TRSSolver::Size,   typename ConvergenceCriterion::Size>::value  );
	
	using Scalar    = typename TRSSolver::Scalar;
	using Size      = typename TRSSolver::Size;
	using Criterion = ConvergenceCriterion;
};
	
template<typename TRSSolver, typename ConvergenceCriterion = L2Norm<typename TRSSolver::Scalar> >
class LSR1TrustRegionSolver 
	: public FirstOrderSolverBase< LSR1TrustRegionSolver<TRSSolver, ConvergenceCriterion> >
	, public TrustRegionMethodBase<typename TRSSolver::Scalar, typename TRSSolver::Size>
	, public LimitedMemorySolverBase<typename TRSSolver::Size>
{
	using Self = LSR1TrustRegionSolver<TRSSolver, ConvergenceCriterion>;
public:
	LNOT_DEFINE_FIRST_ORDER_SOLVER
	LNOT_DEFINE_TRUST_REGION_SOLVER
	LNOT_DEFINE_LIMITED_MEMORY_SOLVER
	
	template<CFirstOrderOracle Oracle, typename ABool> 
	void solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value);
	
	const TRSSolver& getSubproblemSolver() const { return m_trsSolver; }
	      TRSSolver& getSubproblemSolver()       { return m_trsSolver; }
protected:
	LNOT_FIRST_ORDER_SOLVER_ATTRIBUTE
	LNOT_TRUST_REGION_SOLVER_ATTRIBUTE
	LNOT_LIMITED_MEMORY_SOLVER_ATTRIBUTE
private:
	TRSSolver m_trsSolver;
	
	std::unique_ptr<Scalar[]> m_gk;
	std::unique_ptr<Scalar[]> m_gkp1;
	std::unique_ptr<Scalar[]> m_xTrial;
	std::unique_ptr<Scalar[]> m_P;
	std::unique_ptr<Scalar[]> m_Y;
	std::unique_ptr<Scalar[]> m_S;
};

} // namespace LNOT

#include <LNOT/GeneralSolvers/LSR1TrustRegionSolver_impl.hpp>

#endif // LNOT_LSR1_TRUST_REGION_SOLVER_HPP
