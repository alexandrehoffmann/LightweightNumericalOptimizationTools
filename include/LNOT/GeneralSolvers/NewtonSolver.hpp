#ifndef LNOT_NEWTON_SOLVER_HPP
#define LNOT_NEWTON_SOLVER_HPP

#include <LNOT/GeneralSolvers/SecondOrderSolverBase.hpp>
#include <LNOT/LinearSolvers/LinearSolverBase.hpp>
#include <LNOT/LineSearches/LineSearchBase.hpp>

namespace LNOT
{

template<typename LinSolver, typename LineSearch> class NewtonSolver;

template<typename LinSolver, typename LineSearch> 
struct SecondOrderSolverTraits< NewtonSolver<LinSolver, LineSearch> >
{
	static_assert(LinearSolver_concept<LinSolver>);
	static_assert(LineSearch_concept<LineSearch>);
	static_assert(std::is_same<typename LinSolver::Scalar, typename LineSearch::Scalar>::value);
	static_assert(std::is_same<typename LinSolver::Size,   typename LineSearch::Size>::value);
	
	using Scalar = typename LinSolver::Scalar;
	using Size   = typename LinSolver::Size;
};

template<typename LinSolver, typename LineSearch>
class NewtonSolver : public SecondOrderSolverBase< NewtonSolver<LinSolver, LineSearch> >
{
	using Base = SecondOrderSolverBase< NewtonSolver<LinSolver, LineSearch> >;
public:
	using Scalar = typename Base::Scalar;
	using Size   = typename Base::Size;
	using Info   = typename Base::Info;
	
	NewtonSolver(LinSolver& linSolver, LineSearch& lineSearch, const Size maxIt = 200000, const Scalar tol = std::numeric_limits<Scalar>::epsilon()) : Base(maxIt, tol), m_linsSolver(linSolver), m_lineSearch(lineSearch) {}
	
	void clearWorkSpaceImpl();
	
	template<SecondOrderOracle_concept Oracle, bool solveInPlace> 
	void solveImpl(Oracle& oracle, std::bool_constant<solveInPlace> bc, Scalar* x);
private:
	Scalar* m_gk = nullptr;
	Scalar* m_sk = nullptr;
	
	LinSolver&  m_linsSolver;
	LineSearch& m_lineSearch;
};

template<LinearSolver_concept LinSolver, LineSearch_concept LineSearch>
NewtonSolver<LinSolver, LineSearch> makeNewtonSolver(LinSolver& linSolver, LineSearch& lineSearch, const typename LinSolver::Size maxIt = 200000, const typename LinSolver::Scalar tol = std::numeric_limits<typename LinSolver::Scalar>::epsilon()) { return NewtonSolver<LinSolver, LineSearch>(linSolver, lineSearch, maxIt, tol); }

} // namespace LNOT

#include <LNOT/GeneralSolvers/NewtonSolver_impl.hpp>

#endif // LNOT_NEWTON_SOLVER_HPP
