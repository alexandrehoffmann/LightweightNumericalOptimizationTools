#ifndef LNOT_GRADIENT_DESCENT_HPP
#define LNOT_GRADIENT_DESCENT_HPP

#include <LNOT/GeneralSolvers/FirstOrderSolverBase.hpp>
#include <LNOT/LineSearches/LineSearchBase.hpp>
#include <LNOT/ConvergenceCriterions/L2Norm.hpp>

#include <memory>

namespace LNOT
{

template<typename LineSearch, typename ConvergenceCriterion> class GradientDescent;

template<typename LineSearch, typename ConvergenceCriterion> 
struct FirstOrderSolverTraits< GradientDescent<LineSearch, ConvergenceCriterion> >
{
	static_assert(CLineSearch<LineSearch>);
	static_assert(CConvergenceCriterion<ConvergenceCriterion>);
	
	static_assert(std::is_same<typename LineSearch::Scalar, typename ConvergenceCriterion::Scalar>::value);
	static_assert(std::is_same<typename LineSearch::Size,   typename ConvergenceCriterion::Size>::value  );
	
	using Scalar    = typename LineSearch::Scalar;
	using Size      = typename LineSearch::Size;
	using Criterion = ConvergenceCriterion;
};

template<typename LineSearch, typename ConvergenceCriterion = L2Norm<typename LineSearch::Scalar> > 
class GradientDescent : public FirstOrderSolverBase< GradientDescent<LineSearch, ConvergenceCriterion> >
{
	using Self = GradientDescent<LineSearch, ConvergenceCriterion>;
public:
	LNOT_DEFINE_FIRST_ORDER_SOLVER
	
	template<CFirstOrderOracle Oracle, typename ABool> 
	void solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value);
	
	const LineSearch& getLinesearch() const { return m_lineSearch; }
	      LineSearch& getLinesearch()       { return m_lineSearch; }
protected:
	LNOT_FIRST_ORDER_SOLVER_ATTRIBUTE
private:
	std::unique_ptr<Scalar[]> m_gk;
	std::unique_ptr<Scalar[]> m_sk;
	
	LineSearch m_lineSearch;
};

} // namespace LNOT

#include <LNOT/GeneralSolvers/GradientDescent_impl.hpp>

#endif // LNOT_GRADIENT_DESCENT_HPP
