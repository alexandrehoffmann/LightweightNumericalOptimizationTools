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
	
	void clearWorkSpaceImpl();
	
	template<FirstOrderOracle_concept Oracle, bool solveInPlace> 
	void solveImpl(Oracle& oracle, std::bool_constant<solveInPlace> bc, Scalar* x);
	
	const TRSSolver& getSubproblemSolver() const { return m_trsSolver; }
	      TRSSolver& getSubproblemSolver()       { return m_trsSolver; }
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
