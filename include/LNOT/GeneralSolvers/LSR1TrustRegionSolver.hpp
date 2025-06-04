#ifndef LNOT_LSR1_TRUST_REGION_SOLVER_HPP
#define LNOT_LSR1_TRUST_REGION_SOLVER_HPP

#include <LNOT/GeneralSolvers/FirstOrderSolverBase.hpp>
#include <LNOT/GeneralSolvers/TrustRegionMethodBase.hpp>
#include <LNOT/TRSSolvers/TRSSolverBase.hpp>

namespace LNOT
{

template<typename TRSSolver> class LSR1TrustRegionSolver;

template<typename TRSSolver>
struct FirstOrderSolverTraits< LSR1TrustRegionSolver<TRSSolver> >
{
	static_assert(TRSSolver_concept<TRSSolver>);
	
	using Scalar = typename TRSSolver::Scalar;
	using Size   = typename TRSSolver::Size;
};
	
template<typename TRSSolver>
class LSR1TrustRegionSolver 
	: public FirstOrderSolverBase< LSR1TrustRegionSolver<TRSSolver> >
	, public TrustRegionMethodBase<typename TRSSolver::Scalar, typename TRSSolver::Size>
{
	using Base = FirstOrderSolverBase< LSR1TrustRegionSolver<TRSSolver> >;
	using TRSBase = TrustRegionMethodBase<typename TRSSolver::Scalar, typename TRSSolver::Size>;
public:
	using Scalar = typename Base::Scalar;
	using Size   = typename Base::Size;
	using Info   = typename Base::Info;
	
	LSR1TrustRegionSolver(TRSSolver& trsSolver, Size memory, const Size maxIt = 200000, const Scalar tol = std::numeric_limits<Scalar>::epsilon()) : Base(maxIt, tol), m_trsSolver(trsSolver), m_memory(memory) {}
	
	void clearWorkSpace();
	
	template<FirstOrderOracle_concept Oracle, bool solveInPlace> 
	void solve_impl(Oracle& oracle, std::bool_constant<solveInPlace> bc, Scalar* x);
private:
	TRSSolver& m_trsSolver;
	
	Size m_memory;
	
	Scalar* m_gk     = nullptr;
	Scalar* m_gkp1   = nullptr;
	Scalar* m_xTrial = nullptr;
	Scalar* m_P      = nullptr;
	Scalar* m_Y      = nullptr;
	Scalar* m_S      = nullptr;
};

template<TRSSolver_concept TRSSolver>
LSR1TrustRegionSolver<TRSSolver> makeLSR1Solver(TRSSolver& trsSolver, const typename TRSSolver::Size memory, const typename TRSSolver::Size maxIt = 200000, const typename TRSSolver::Scalar tol = std::numeric_limits<typename TRSSolver::Scalar>::epsilon()) { return LSR1TrustRegionSolver<TRSSolver>(trsSolver, memory, maxIt, tol); }

} // namespace LNOT

#include <LNOT/GeneralSolvers/LSR1TrustRegionSolver_impl.hpp>

#endif // LNOT_LSR1_TRUST_REGION_SOLVER_HPP
