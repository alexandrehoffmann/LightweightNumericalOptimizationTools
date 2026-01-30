#ifndef LNOT_BFGS_HPP
#define LNOT_BFGS_HPP

#include <LNOT/GeneralSolvers/FirstOrderSolverBase.hpp>
#include <LNOT/LineSearches/LineSearchBase.hpp>

namespace LNOT
{

template<typename LineSearch> class BFGS;

template<typename LineSearch> 
struct FirstOrderSolverTraits< BFGS<LineSearch> >
{
	static_assert(CLineSearch<LineSearch>);
	
	using Scalar = typename LineSearch::Scalar;
	using Size   = typename LineSearch::Size;
};

template<typename LineSearch>
class BFGS : public FirstOrderSolverBase< BFGS<LineSearch> >
{
	using Self = BFGS<LineSearch>;
public:
	LNOT_DEFINE_FIRST_ORDER_SOLVER
	
	void clearWorkSpaceImpl();
	
	template<CFirstOrderOracle Oracle, typename ABool> 
	void solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value);
	
	const LineSearch& getLinesearch() const { return m_lineSearch; }
	      LineSearch& getLinesearch()       { return m_lineSearch; }
protected:
	LNOT_FIRST_ORDER_SOLVER_ATTRIBUTE
private:
	Scalar* m_gk    = nullptr;
	Scalar* m_gkp1  = nullptr;
	Scalar* m_sk    = nullptr;
	Scalar* m_yk    = nullptr;
	Scalar* m_uk    = nullptr;
	Scalar* m_invBk = nullptr;
	
	LineSearch m_lineSearch;
};

} // namespace LNOT

#include <LNOT/GeneralSolvers/BFGS_impl.hpp>

#endif // LNOT_BFGS_HPP
