#ifndef LNOT_NEWTON_TRUST_REGION_SOLVER_HPP
#define LNOT_NEWTON_TRUST_REGION_SOLVER_HPP

#include <LNOT/GeneralSolvers/SecondOrderSolverBase.hpp>
#include <LNOT/GeneralSolvers/TrustRegionMethodBase.hpp>
#include <LNOT/TRSSolvers/TRSSolverBase.hpp>
#include <LNOT/ConvergenceCriterions/L2Norm.hpp>
#include <LNOT/Traits.hpp>

#include <memory>

namespace LNOT
{
	
template<typename TRSSolver, typename ConvergenceCriterion> class NewtonTrustRegionSolver;

template<typename TRSSolver, typename ConvergenceCriterion>
struct SecondOrderSolverTraits< NewtonTrustRegionSolver<TRSSolver, ConvergenceCriterion> >
{
	static_assert(CTRSSolver<TRSSolver>);
	static_assert(CConvergenceCriterion<ConvergenceCriterion>);
	
	static_assert(std::same_as< ScalarFor<TRSSolver>, ScalarFor<ConvergenceCriterion> >);
	static_assert(std::same_as< SizeFor<TRSSolver>,   SizeFor<ConvergenceCriterion> >);
	
	using Scalar    = ScalarFor<TRSSolver>;
	using Size      = SizeFor<TRSSolver>;
	using Criterion = ConvergenceCriterion;
};
	
template<typename TRSSolver, typename ConvergenceCriterion = L2Norm<ScalarFor<TRSSolver>>>
class NewtonTrustRegionSolver 
	: public SecondOrderSolverBase< NewtonTrustRegionSolver<TRSSolver, ConvergenceCriterion> >
	, public TrustRegionMethodBase<ScalarFor<TRSSolver>, SizeFor<TRSSolver>>
{
	using Self = NewtonTrustRegionSolver<TRSSolver, ConvergenceCriterion>;
public:
	LNOT_DEFINE_SECOND_ORDER_SOLVER
	LNOT_DEFINE_TRUST_REGION_SOLVER
	
	template<CSecondOrderOracle Oracle, typename ABool> 
	void solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value);
	
	const TRSSolver& getSubproblemSolver() const { return m_trsSolver; }
	      TRSSolver& getSubproblemSolver()       { return m_trsSolver; }
protected:
	LNOT_SECOND_ORDER_SOLVER_ATTRIBUTE
	LNOT_TRUST_REGION_SOLVER_ATTRIBUTE
private:
	TRSSolver m_trsSolver;
	
	std::unique_ptr<Scalar[]> m_gk;
	std::unique_ptr<Scalar[]> m_sk;
	std::unique_ptr<Scalar[]> m_xTrial;
};

} // namespace LNOT

#include <LNOT/GeneralSolvers/NewtonTrustRegionSolver_impl.hpp>

#endif // LNOT_NEWTON_TRUST_REGION_SOLVER_HPP
