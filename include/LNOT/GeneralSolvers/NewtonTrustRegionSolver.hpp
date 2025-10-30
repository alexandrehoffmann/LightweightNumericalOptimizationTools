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
	static_assert(TRSSolver_concept<TRSSolver>);
	
	using Scalar = typename TRSSolver::Scalar;
	using Size   = typename TRSSolver::Size;
};
	
template<typename TRSSolver>
class NewtonTrustRegionSolver 
	: public SecondOrderSolverBase< NewtonTrustRegionSolver<TRSSolver> >
	, public TrustRegionMethodBase<typename TRSSolver::Scalar, typename TRSSolver::Size>
{
	using Base = SecondOrderSolverBase< NewtonTrustRegionSolver<TRSSolver> >;
	using TRSBase = TrustRegionMethodBase<typename TRSSolver::Scalar, typename TRSSolver::Size>;
public:
	using Scalar = typename Base::Scalar;
	using Size   = typename Base::Size;
	using Info   = typename Base::Info;
	
	NewtonTrustRegionSolver(TRSSolver& trsSolver, const Size maxIt = 200000, const Scalar tol = std::numeric_limits<Scalar>::epsilon()) : Base(maxIt, tol), m_trsSolver(trsSolver) {}
	
	void clearWorkSpaceImpl();
	
	template<SecondOrderOracle_concept Oracle, bool solveInPlace> 
	void solveImpl(Oracle& oracle, std::bool_constant<solveInPlace> bc, Scalar* x);
private:
	TRSSolver& m_trsSolver;
	
	Scalar* m_gk     = nullptr;
	Scalar* m_sk     = nullptr;
	Scalar* m_xTrial = nullptr;
};

template<TRSSolver_concept TRSSolver>
NewtonTrustRegionSolver<TRSSolver> makeNewtonSolver(TRSSolver& trsSolver, const typename TRSSolver::Size maxIt = 200000, const typename TRSSolver::Scalar tol = std::numeric_limits<typename TRSSolver::Scalar>::epsilon()) { return NewtonTrustRegionSolver<TRSSolver>(trsSolver, maxIt, tol); }

} // namespace LNOT

#include <LNOT/GeneralSolvers/NewtonTrustRegionSolver_impl.hpp>

#endif // LNOT_NEWTON_TRUST_REGION_SOLVER_HPP
