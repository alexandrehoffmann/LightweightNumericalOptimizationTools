#ifndef LNOT_LBFGS_HPP
#define LNOT_LBFGS_HPP

#include <LNOT/GeneralSolvers/FirstOrderSolverBase.hpp>
#include <LNOT/GeneralSolvers/LimitedMemorySolverBase.hpp>
#include <LNOT/LineSearches/LineSearchBase.hpp>
#include <LNOT/ConvergenceCriterions/L2Norm.hpp>

#include <memory>

namespace LNOT
{

template<typename LineSearch, typename ConvergenceCriterion> class LBFGS;

template<typename LineSearch, typename ConvergenceCriterion> 
struct FirstOrderSolverTraits< LBFGS<LineSearch, ConvergenceCriterion> >
{
	static_assert(CLineSearch<LineSearch>);
	static_assert(CConvergenceCriterion<ConvergenceCriterion>);
	
	static_assert(std::same_as< LineSearch_Scalar<LineSearch>, Criterion_Scalar<ConvergenceCriterion> >);
	static_assert(std::same_as< LineSearch_Size<LineSearch>,   Criterion_Size<ConvergenceCriterion> >);
	
	using Scalar    = LineSearch_Scalar<LineSearch>;
	using Size      = LineSearch_Size<LineSearch>;
	using Criterion = ConvergenceCriterion;
};

template<typename LineSearch, typename ConvergenceCriterion = L2Norm<LineSearch_Scalar<LineSearch>>>
class LBFGS 
	: public FirstOrderSolverBase< LBFGS<LineSearch, ConvergenceCriterion> >
	, public LimitedMemorySolverBase< LineSearch_Size<LineSearch> >
{
	using Self = LBFGS<LineSearch, ConvergenceCriterion>;
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
