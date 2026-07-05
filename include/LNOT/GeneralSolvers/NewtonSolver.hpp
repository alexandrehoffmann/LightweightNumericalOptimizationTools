#ifndef LNOT_NEWTON_SOLVER_HPP
#define LNOT_NEWTON_SOLVER_HPP

#include <LNOT/GeneralSolvers/SecondOrderSolverBase.hpp>
#include <LNOT/LinearSolvers/LinearSolverBase.hpp>
#include <LNOT/LineSearches/LineSearchBase.hpp>
#include <LNOT/ConvergenceCriterions/L2Norm.hpp>
#include <LNOT/Traits.hpp>

#include <memory>

namespace LNOT
{

template<typename LinSolver, typename LineSearch, typename ConvergenceCriterion> class NewtonSolver;

template<typename LinSolver, typename LineSearch, typename ConvergenceCriterion> 
struct SecondOrderSolverTraits< NewtonSolver<LinSolver, LineSearch, ConvergenceCriterion> >
{
	static_assert(CLinearSolver<LinSolver>);
	static_assert(CLineSearch<LineSearch>);
	
	static_assert(CConvergenceCriterion<ConvergenceCriterion>);
	static_assert(std::same_as< ScalarFor<LinSolver>, ScalarFor<LineSearch> > and std::same_as< ScalarFor<LinSolver>, ScalarFor<ConvergenceCriterion> >);
	static_assert(std::same_as< SizeFor<LinSolver>,   SizeFor<LineSearch> >   and std::same_as< SizeFor<LinSolver>,   SizeFor<ConvergenceCriterion> >  );
	
	using Scalar    = ScalarFor<LinSolver>;
	using Size      = SizeFor<LinSolver>;
	using Criterion = ConvergenceCriterion;
};

template<typename LinSolver, typename LineSearch, typename ConvergenceCriterion = L2Norm<ScalarFor<LinSolver>>>
class NewtonSolver : public SecondOrderSolverBase< NewtonSolver<LinSolver, LineSearch, ConvergenceCriterion> >
{
	using Self = NewtonSolver<LinSolver, LineSearch, ConvergenceCriterion>;
public:
	LNOT_DEFINE_SECOND_ORDER_SOLVER
	
	template<CSecondOrderOracle Oracle, typename ABool> 
	void solveImpl(Oracle&& oracle, const ABool solveInPlace, Scalar* x) requires(isBool<ABool>);
	
	const LinSolver& getLinearSolver() const { return m_linsSolver; }
	      LinSolver& getLinearSolver()       { return m_linsSolver; }
	
	const LineSearch& getLinesearch() const { return m_lineSearch; }
	      LineSearch& getLinesearch()       { return m_lineSearch; }
protected:
	LNOT_SECOND_ORDER_SOLVER_ATTRIBUTE
private:
	std::unique_ptr<Scalar[]> m_gk;
	std::unique_ptr<Scalar[]> m_sk;
	
	LinSolver  m_linsSolver;
	LineSearch m_lineSearch;
};

} // namespace LNOT

#include <LNOT/GeneralSolvers/NewtonSolver_impl.hpp>

#endif // LNOT_NEWTON_SOLVER_HPP
