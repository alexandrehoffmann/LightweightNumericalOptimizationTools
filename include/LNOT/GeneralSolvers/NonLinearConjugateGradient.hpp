#ifndef LNOT_NON_LINEAR_CONJUGATE_GRADIENT_HPP
#define LNOT_NON_LINEAR_CONJUGATE_GRADIENT_HPP

#include <LNOT/GeneralSolvers/FirstOrderSolverBase.hpp>
#include <LNOT/GeneralSolvers/NLCGUpdateStrategy.hpp>
#include <LNOT/LineSearches/LineSearchBase.hpp>
#include <LNOT/ConvergenceCriterions/L2Norm.hpp>

#include <memory>

namespace LNOT
{

template<typename LineSearch, NLCGUpdateStrategy UpdateStrategy, typename ConvergenceCriterion> class NonLinearConjugateGradient;

template<typename LineSearch, NLCGUpdateStrategy UpdateStrategy, typename ConvergenceCriterion> 
struct FirstOrderSolverTraits< NonLinearConjugateGradient<LineSearch, UpdateStrategy, ConvergenceCriterion> >
{
	static_assert(CLineSearch<LineSearch>);
	static_assert(CConvergenceCriterion<ConvergenceCriterion>);
	
	static_assert(std::same_as< LineSearch_Scalar<LineSearch>, Criterion_Scalar<ConvergenceCriterion> >);
	static_assert(std::same_as< LineSearch_Size<LineSearch>,   Criterion_Size<ConvergenceCriterion> >);
	
	using Scalar    = LineSearch_Scalar<LineSearch>;
	using Size      = LineSearch_Size<LineSearch>;
	using Criterion = ConvergenceCriterion;
};

template<typename LineSearch, NLCGUpdateStrategy UpdateStrategy = NLCGUpdateStrategy::HESTENES_STIEFEL, typename ConvergenceCriterion = L2Norm<LineSearch_Scalar<LineSearch>>> 
class NonLinearConjugateGradient : public FirstOrderSolverBase< NonLinearConjugateGradient<LineSearch, UpdateStrategy, ConvergenceCriterion> >
{
	using Self = NonLinearConjugateGradient<LineSearch, UpdateStrategy, ConvergenceCriterion>;
public:
	LNOT_DEFINE_FIRST_ORDER_SOLVER
	
	template<CFirstOrderOracle Oracle, typename ABool> 
	void solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value);
	
	Scalar getBeta(const Size size);
	
	const LineSearch& getLinesearch() const { return m_lineSearch; }
	      LineSearch& getLinesearch()       { return m_lineSearch; }
protected:
	LNOT_FIRST_ORDER_SOLVER_ATTRIBUTE
private:
	std::unique_ptr<Scalar[]> m_dk;
	std::unique_ptr<Scalar[]> m_yk;
	std::unique_ptr<Scalar[]> m_gk;
	std::unique_ptr<Scalar[]> m_gkp1;
	
	LineSearch m_lineSearch;
};

} // namespace LNOT

#include <LNOT/GeneralSolvers/NonLinearConjugateGradient_impl.hpp>

#endif // LNOT_NON_LINEAR_CONJUGATE_GRADIENT_HPP
