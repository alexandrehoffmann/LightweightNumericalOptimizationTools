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
	using Base = FirstOrderSolverBase< SR1TrustRegionSolver<TRSSolver> >;
	using TRSBase = TrustRegionMethodBase<typename TRSSolver::Scalar, typename TRSSolver::Size>;
public:
	using Scalar = typename Base::Scalar;
	using Size   = typename Base::Size;
	using Info   = typename Base::Info;
	
	SR1TrustRegionSolver(TRSSolver& trsSolver, const Size maxIt = 200000, const Scalar tol = std::numeric_limits<Scalar>::epsilon()) : Base(maxIt, tol), m_trsSolver(trsSolver) {}
	
	void clearWorkSpace();
	
	template<FirstOrderOracle_concept Oracle, bool solveInPlace> 
	void solve_impl(Oracle& oracle, std::bool_constant<solveInPlace> bc, Scalar* x);
private:
	TRSSolver& m_trsSolver;
	
	Scalar* m_gk     = nullptr;
	Scalar* m_gkp1   = nullptr;
	Scalar* m_uk     = nullptr;
	Scalar* m_sk     = nullptr;
	Scalar* m_Bk     = nullptr;
	Scalar* m_xTrial = nullptr;
};

template<TRSSolver_concept TRSSolver>
SR1TrustRegionSolver<TRSSolver> makeSR1Solver(TRSSolver& trsSolver, const typename TRSSolver::Size maxIt = 200000, const typename TRSSolver::Scalar tol = std::numeric_limits<typename TRSSolver::Scalar>::epsilon()) { return SR1TrustRegionSolver<TRSSolver>(trsSolver, maxIt, tol); }

} // namespace LNOT

#include <LNOT/GeneralSolvers/SR1TrustRegionSolver_impl.hpp>

#endif // LNOT_SR1_TRUST_REGION_SOLVER_HPP
