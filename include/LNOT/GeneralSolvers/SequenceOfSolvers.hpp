#ifndef LNOT_SEQUENCE_OF_SOLVERS_HPP
#define LNOT_SEQUENCE_OF_SOLVERS_HPP

#include <LNOT/GeneralSolvers/FirstOrderSolverBase.hpp>
#include <LNOT/GeneralSolvers/SecondOrderSolverBase.hpp>

#include <BIC/Core.hpp>

#include <array>
#include <tuple>

namespace LNOT
{
	
template<class Solver>
concept CSolver = CFirstOrderSolver<Solver> or CSecondOrderSolver<Solver>;

template<CSolver... Solvers> class SequenceOfSolvers;

template<CSolver... Solvers> requires(std::max({Solvers::order...}) == 1)
struct FirstOrderSolverTraits< SequenceOfSolvers<Solvers...> >
{
	using FirstSolver = std::tuple_element_t<0, std::tuple<Solvers...>>;
	
	static_assert((std::is_same<typename FirstSolver::Scalar, typename Solvers::Scalar>::value and ...));
	static_assert((std::is_same<typename FirstSolver::Size,   typename Solvers::Size>::value and ...));
	
	using Scalar = typename FirstSolver::Scalar;
	using Size   = typename FirstSolver::Size;
};

template<CSolver... Solvers> requires(std::max({Solvers::order...}) == 2)
struct SecondOrderSolverTraits< SequenceOfSolvers<Solvers...> >
{
	using FirstSolver = std::tuple_element_t<0, std::tuple<Solvers...>>;
	
	static_assert((std::is_same<typename FirstSolver::Scalar, typename Solvers::Scalar>::value and ...));
	static_assert((std::is_same<typename FirstSolver::Size,   typename Solvers::Size>::value and ...));
	
	using Scalar = typename FirstSolver::Scalar;
	using Size   = typename FirstSolver::Size;
};

namespace detail
{

template<CSolver... Solvers>
using SequenceOfSolversBase = std::conditional_t<std::max({Solvers::order...}) == 2
	, SecondOrderSolverBase< SequenceOfSolvers<Solvers...> >
	, FirstOrderSolverBase< SequenceOfSolvers<Solvers...> > 
>;

} // namespace detail

template<CSolver... Solvers>
class SequenceOfSolvers : public detail::SequenceOfSolversBase<Solvers...>
{
	static constexpr inline size_t nSolvers = sizeof...(Solvers);
public:
	using Base   = detail::SequenceOfSolversBase<Solvers...>;
	using Size   = typename Base::Size;
	using Scalar = typename Base::Scalar;
	using Info   = typename Base::Info;
	
	template<typename ABool> using IsBool = typename Base::template IsBool<ABool>; 
	template<typename ASize> using IsSize = typename Base::template IsSize<ASize>; 
	
	template<size_t I> using IthSolver = std::tuple_element_t<I, std::tuple<Solvers...>>;
	
	SequenceOfSolvers(const Size maxIt = 200000, const Scalar relTol = Base::defaultEps, const Scalar absTol = Base::defaultEps);
	
	SequenceOfSolvers(const std::array<Size, nSolvers>& maxIt, const std::array<Scalar, nSolvers>& relTol, const std::array<Scalar, nSolvers>& absTol);
	
	SequenceOfSolvers(Solvers&&... solver) : m_solvers(std::move(solver)...) {}
	
	void setMaxIt  (const std::array<Size  , nSolvers>& maxIt);
	void setRelTol (const std::array<Scalar, nSolvers>& tol);
	void setAbsTol (const std::array<Scalar, nSolvers>& tol);
	void setTol    (const std::array<Scalar, nSolvers>& tol) { setRelTol(tol); setAbsTol(tol); }
	
	template<CFirstOrderOracle Oracle, typename ABool> 
	void solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value and Base::order == 1) { solveImplDetail(oracle, solveInPlace, x); }
	
	template<CSecondOrderOracle Oracle, typename ABool> 
	void solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value and Base::order == 2) { solveImplDetail(oracle, solveInPlace, x); }
	
	template<size_t I> const IthSolver<I>& getSolver() const { return std::get<I>(m_solvers); }
	template<size_t I>       IthSolver<I>& getSolver()       { return std::get<I>(m_solvers); }
	
	template<size_t I> const IthSolver<I>& getSolver(const BIC::Fixed<size_t, I> i) const { return getSolver<i>(); }
	template<size_t I>       IthSolver<I>& getSolver(const BIC::Fixed<size_t, I> i)       { return getSolver<i>(); }
protected:
	using Base::m_maxIt;
	using Base::m_relTol;
	using Base::m_absTol;
	using Base::m_nIt;
	using Base::m_info;
	using Base::m_fx;
	using Base::m_residual;
	using Base::m_workCapacity;
	using Base::m_innerIts;
	using Base::m_out;
private:
	template<class Oracle, typename ABool> 
	void solveImplDetail(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value);

	std::tuple<Solvers...> m_solvers;
};
	
} // namespace LNOT

#include <LNOT/GeneralSolvers/SequenceOfSolvers_impl.hpp>

#endif // LNOT_SEQUENCE_OF_SOLVERS_HPP
