#ifndef LNOT_LBFGS_HPP
#define LNOT_LBFGS_HPP

#include <LNOT/GeneralSolvers/FirstOrderSolverBase.hpp>
#include <LNOT/GeneralSolvers/LimitedMemorySolverBase.hpp>
#include <LNOT/LineSearches/LineSearchBase.hpp>

#include <memory>

namespace LNOT
{

template<typename LineSearch> class LBFGS;

template<typename LineSearch> 
struct FirstOrderSolverTraits< LBFGS<LineSearch> >
{
	static_assert(CLineSearch<LineSearch>);
	
	using Scalar = typename LineSearch::Scalar;
	using Size   = typename LineSearch::Size;
};

template<typename LineSearch>
class LBFGS 
	: public FirstOrderSolverBase< LBFGS<LineSearch> >
	, public LimitedMemorySolverBase<typename LineSearch::Size>
{
	using Self = LBFGS<LineSearch>;
public:
	LNOT_DEFINE_FIRST_ORDER_SOLVER
	LNOT_DEFINE_LIMITED_MEMORY_SOLVER
	
	template<CFirstOrderOracle Oracle, typename ABool> 
	void solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value);
	
	const LineSearch& getLinesearch() const { return m_lineSearch; }
	      LineSearch& getLinesearch()       { return m_lineSearch; }
protected:
	LNOT_FIRST_ORDER_SOLVER_ATTRIBUTE
	LNOT_LIMITED_MEMORY_SOLVER_ATTRIBUTE
private:
	std::unique_ptr<Scalar[]> m_gk;
	std::unique_ptr<Scalar[]> m_gkp1;
	std::unique_ptr<Scalar[]> m_S;
	std::unique_ptr<Scalar[]> m_Y;
	std::unique_ptr<Scalar[]> m_dk;
	
	LineSearch m_lineSearch;
};

} // namespace LNOT

#include <LNOT/GeneralSolvers/LBFGS_impl.hpp>

#endif // LNOT_LBFGS_HPP
