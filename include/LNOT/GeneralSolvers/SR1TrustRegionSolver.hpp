#ifndef LNOT_SR1_TRUST_REGION_SOLVER_HPP
#define LNOT_SR1_TRUST_REGION_SOLVER_HPP

#include <LNOT/GeneralSolvers/FirstOrderSolverBase.hpp>
#include <LNOT/GeneralSolvers/TrustRegionMethodBase.hpp>
#include <LNOT/TRSSolvers/TRSSolverBase.hpp>
#include <LNOT/ConvergenceCriterions/L2Norm.hpp>

#include <memory>

namespace LNOT
{

template<typename TRSSolver, typename ConvergenceCriterion> class SR1TrustRegionSolver;

template<typename TRSSolver, typename ConvergenceCriterion>
struct FirstOrderSolverTraits< SR1TrustRegionSolver<TRSSolver, ConvergenceCriterion> >
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
class SR1TrustRegionSolver 
	: public FirstOrderSolverBase< SR1TrustRegionSolver<TRSSolver, ConvergenceCriterion> >
	, public TrustRegionMethodBase<typename TRSSolver::Scalar, typename TRSSolver::Size>
{
	using Self = SR1TrustRegionSolver<TRSSolver, ConvergenceCriterion>;
	using TRM  = TrustRegionMethodBase<typename TRSSolver::Scalar, typename TRSSolver::Size>;
public:
	LNOT_DEFINE_FIRST_ORDER_SOLVER
	LNOT_DEFINE_TRUST_REGION_SOLVER
	
	template<CFirstOrderOracle Oracle, typename ABool> 
	void solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value);
	
	const TRSSolver& getSubproblemSolver() const { return m_trsSolver; }
	      TRSSolver& getSubproblemSolver()       { return m_trsSolver; }
protected:
	LNOT_FIRST_ORDER_SOLVER_ATTRIBUTE
	LNOT_TRUST_REGION_SOLVER_ATTRIBUTE
private:
	TRSSolver m_trsSolver;
	
	std::unique_ptr<Scalar[]> m_gk;
	std::unique_ptr<Scalar[]> m_gkp1;
	std::unique_ptr<Scalar[]> m_uk;
	std::unique_ptr<Scalar[]> m_sk;
	std::unique_ptr<Scalar[]> m_Bk;
	std::unique_ptr<Scalar[]> m_xTrial;
};

} // namespace LNOT

#include <LNOT/GeneralSolvers/SR1TrustRegionSolver_impl.hpp>

#endif // LNOT_SR1_TRUST_REGION_SOLVER_HPP
