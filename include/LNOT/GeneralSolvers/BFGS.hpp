#ifndef LNOT_BFGS_HPP
#define LNOT_BFGS_HPP

#include <LNOT/GeneralSolvers/FirstOrderSolverBase.hpp>
#include <LNOT/LineSearches/LineSearchBase.hpp>
#include <LNOT/ConvergenceCriterions/L2Norm.hpp>

#include <memory>

namespace LNOT
{

template<typename LineSearch, typename ConvergenceCriterion> class BFGS;

template<typename LineSearch, typename ConvergenceCriterion> 
struct FirstOrderSolverTraits< BFGS<LineSearch, ConvergenceCriterion> >
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
class BFGS : public FirstOrderSolverBase< BFGS<LineSearch, ConvergenceCriterion> >
{
	using Self = BFGS<LineSearch, ConvergenceCriterion>;
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
	std::unique_ptr<Scalar[]> m_gkp1;
	std::unique_ptr<Scalar[]> m_sk;
	std::unique_ptr<Scalar[]> m_yk;
	std::unique_ptr<Scalar[]> m_uk;
	std::unique_ptr<Scalar[]> m_invBk;
	
	LineSearch m_lineSearch;
};

} // namespace LNOT

#include <LNOT/GeneralSolvers/BFGS_impl.hpp>

#endif // LNOT_BFGS_HPP
