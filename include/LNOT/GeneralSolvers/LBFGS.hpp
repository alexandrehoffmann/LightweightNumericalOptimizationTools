#ifndef LNOT_LBFGS_HPP
#define LNOT_LBFGS_HPP

#include <LNOT/GeneralSolvers/FirstOrderSolverBase.hpp>
#include <LNOT/GeneralSolvers/LimitedMemorySolverBase.hpp>
#include <LNOT/LineSearches/LineSearchBase.hpp>

namespace LNOT
{

template<typename LineSearch> class LBFGS;

template<typename LineSearch> 
struct FirstOrderSolverTraits< LBFGS<LineSearch> >
{
	static_assert(LineSearch_concept<LineSearch>);
	
	using Scalar = typename LineSearch::Scalar;
	using Size   = typename LineSearch::Size;
};

template<typename LineSearch>
class LBFGS 
	: public FirstOrderSolverBase< LBFGS<LineSearch> >
	, public LimitedMemorySolverBase<typename LineSearch::Size>
{
	using Self = LBFGS<LineSearch>;
	using LM   = LimitedMemorySolverBase<typename LineSearch::Size>;
public:
	LNOT_DEFINE_FIRST_ORDER_SOLVER
	LNOT_DEFINE_LIMITED_MEMORY_SOLVER
	
	void clearWorkSpaceImpl();
	
	template<FirstOrderOracle_concept Oracle, typename ABool> 
	void solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value);
	
	const LineSearch& getLinesearch() const { return m_lineSearch; }
	      LineSearch& getLinesearch()       { return m_lineSearch; }
protected:
	LNOT_FIRST_ORDER_SOLVER_ATTRIBUTE
	LNOT_LIMITED_MEMORY_SOLVER_ATTRIBUTE
private:
	Scalar* m_gk   = nullptr;
	Scalar* m_gkp1 = nullptr;
	Scalar* m_S    = nullptr;
	Scalar* m_Y    = nullptr;
	Scalar* m_dk   = nullptr;
	
	LineSearch m_lineSearch;
};

} // namespace LNOT

#include <LNOT/GeneralSolvers/LBFGS_impl.hpp>

#endif // LNOT_LBFGS_HPP
