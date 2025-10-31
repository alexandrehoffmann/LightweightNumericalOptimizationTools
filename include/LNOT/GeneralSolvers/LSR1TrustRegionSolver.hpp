#ifndef LNOT_LSR1_TRUST_REGION_SOLVER_HPP
#define LNOT_LSR1_TRUST_REGION_SOLVER_HPP

#include <LNOT/GeneralSolvers/FirstOrderSolverBase.hpp>
#include <LNOT/GeneralSolvers/TrustRegionMethodBase.hpp>
#include <LNOT/GeneralSolvers/LimitedMemorySolverBase.hpp>
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
	, public LimitedMemorySolverBase<typename TRSSolver::Size>
{
	using Base    = FirstOrderSolverBase< LSR1TrustRegionSolver<TRSSolver> >;
	using TRSBase = TrustRegionMethodBase<typename TRSSolver::Scalar, typename TRSSolver::Size>;
	using LMBase  = LimitedMemorySolverBase<typename TRSSolver::Size>;
public:
	using Scalar = typename Base::Scalar;
	using Size   = typename Base::Size;
	using Info   = typename Base::Info;
	
	void clearWorkSpaceImpl();
	
	template<FirstOrderOracle_concept Oracle, bool solveInPlace> 
	void solveImpl(Oracle& oracle, std::bool_constant<solveInPlace> bc, Scalar* x);
	
	const TRSSolver& getSubproblemSolver() const { return m_trsSolver; }
	      TRSSolver& getSubproblemSolver()       { return m_trsSolver; }
private:
	TRSSolver m_trsSolver;
	
	Scalar* m_gk     = nullptr;
	Scalar* m_gkp1   = nullptr;
	Scalar* m_xTrial = nullptr;
	Scalar* m_P      = nullptr;
	Scalar* m_Y      = nullptr;
	Scalar* m_S      = nullptr;
};

} // namespace LNOT

#include <LNOT/GeneralSolvers/LSR1TrustRegionSolver_impl.hpp>

#endif // LNOT_LSR1_TRUST_REGION_SOLVER_HPP
