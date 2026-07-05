#ifndef LNOT_FIRST_ORDER_SOLVER_BASE_HPP
#define LNOT_FIRST_ORDER_SOLVER_BASE_HPP

#include <LNOT/Oracles/OracleWrapper.hpp>
#include <LNOT/Oracles/concepts.hpp>
#include <LNOT/FloatingPoint/NumTraits.hpp>
#include <LNOT/CRTPBase.hpp>

#include <BIC/Core.hpp>

#include <cstdio>
#include <cmath>

#include <vector>
#include <span>
#include <ranges>
#include <algorithm>

namespace LNOT
{

template<class Derived> struct FirstOrderSolverTraits;

/**
 * @brief Base class template for first-order (gradient-based) optimization solvers.
 * 
 * Uses CRTP (Curiously Recurring Template Pattern) to enable static polymorphism for derived solvers.
 * Provides utilities for managing iterations, stopping criteria, and evaluation metrics.
 *
 * @tparam Derived The actual solver class inheriting from this base.
 */
template<class Derived>
class FirstOrderSolverBase : public CRTPBase<Derived>
{
	using DerivedTraits = FirstOrderSolverTraits<Derived>;
	using CRTP          = CRTPBase<Derived>;
public:
	using Scalar    = typename DerivedTraits::Scalar;    ///<  @brief The scalar type used in computations (e.g., float, double)
	using Size      = typename DerivedTraits::Size;      ///<  @brief The size type used for indexing and loop counters
	using Criterion = typename DerivedTraits::Criterion; ///<  @brief The convergence criterion
	
	/// @brief Enumeration indicating solver termination status.
	enum class Info 
	{
		SUCCESS,  ///<  Converged successfully
		FAILURE,  ///<  Failed to converge
		BREAKDOWN ///<  Numerical breakdown (e.g., division by zero)
	};
	
	template<typename ASize> static constexpr bool isSize = std::same_as<Size, BIC::Mutable<ASize>>;
	template<typename ABool> static constexpr bool isBool = std::same_as<bool, BIC::Mutable<ABool>>;
	
	static inline constexpr Size order = 1;
	
	static inline const Scalar defaultEps = NumTraits<Scalar>::epsilon; ///<  @brief Default value for relative and absolute tolerance of the solver.  
	
	// ========================================================================
	// SOLVER INTERFACES
	// ========================================================================
	
	/**
	 * @brief Constructor
	 * @param maxIt Maximum number of iterations (default: 200000).
	 * @param tol Convergence tolerance (default: machine epsilon).
	 */
	FirstOrderSolverBase(const Size maxIt = 200000, const Scalar relTol = defaultEps, const Scalar absTol = defaultEps) : m_maxIt(maxIt), m_relTol(relTol), m_absTol(absTol) {}
	
	/**
	 * @brief Solve using a valid FirstOrderOracle (no initial guess).
	 * @param oracle An oracle object conforming to CFirstOrderOracle.
	 * @param x Output vector for the solution.
	 */
	template<CFirstOrderOracle Oracle>
	void solve(Oracle&& oracle, Scalar* x) { CRTP::derived().solveImpl(std::forward<Oracle>(oracle), BIC::fixed<bool,false>, x); }

	/**
	 * @brief Solve using a valid FirstOrderOracle (no initial guess).
	 * @param oracle An oracle object conforming to CFirstOrderOracle.
	 * @param x Output vector for the solution.
	 */
	template<CFirstOrderOracle Oracle>
	void solve(Oracle&& oracle, const std::span<Scalar> x) { assert(x.size() == oracle.getNDims()); solve(std::forward<Oracle>(oracle), x.data()); }
	
	/**
	 * @brief Solve using a valid FirstOrderOracle in place (uses x as an initial guess).
	 * @param oracle An oracle object conforming to CFirstOrderOracle.
	 * @param x Initial guess and output vector for the solution.
	 */
	template<CFirstOrderOracle Oracle>
	void solveInPlace(Oracle&& oracle, Scalar* x) { CRTP::derived().solveImpl(std::forward<Oracle>(oracle), BIC::fixed<bool,true>, x); }

	/**
	 * @brief Solve using a valid FirstOrderOracle in place (uses x as an initial guess).
	 * @param oracle An oracle object conforming to CFirstOrderOracle.
	 * @param x Initial guess and output vector for the solution.
	 */
	template<CFirstOrderOracle Oracle>
	void solveInPlace(Oracle&& oracle, const std::span<Scalar> x) { assert(x.size() == oracle.getNDims()); solveInPlace(std::forward<Oracle>(oracle), x.data()); }
	
	/**
	 * @brief Solve using a valid FirstOrderOracle with an initial guess.
	 * @param oracle An oracle object.
	 * @param x0 Initial guess.
	 * @param x Output vector for the solution.
	 */
	template<CFirstOrderOracle Oracle>
	void solveWithGuess(Oracle&& oracle, const Scalar* x0, Scalar* x) { std::copy(x0, x0 + oracle.getNDims(), x); solveInPlace(std::forward<Oracle>(oracle), x); }

	/**
	 * @brief Solve using a valid FirstOrderOracle with an initial guess.
	 * @param oracle An oracle object.
	 * @param x0 Initial guess.
	 * @param x Output vector for the solution.
	 */
	template<CFirstOrderOracle Oracle, std::ranges::input_range InitialGuess>
	void solveWithGuess(Oracle&& oracle, InitialGuess&& x0, const std::span<Scalar> x) { std::ranges::copy(x0, x.begin()); solveInPlace(std::forward<Oracle>(oracle), x); }
	
	/**
	 * @brief Solve using raw function and gradient functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param size Problem dimension.
	 * @param x Output solution vector.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient, typename ASize>
	void solve(Function&& f, Gradient&& g, const ASize size, Scalar* x) requires (isSize<ASize>) { solve(makeOracle<Scalar>(size, std::forward<Function>(f), std::forward<Gradient>(g)), x); }

	/**
	 * @brief Solve using raw function and gradient functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param x Output solution vector.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient>
	void solve(Function&& f, Gradient&& g, const std::span<Scalar> x) { solve(makeOracle<Scalar>(Size(x.size()), std::forward<Function>(f), std::forward<Gradient>(g)), x.data()); }
	
	/**
	 * @brief Solve in place using raw function and gradient functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param size Problem dimension.
	 * @param x Initial guess and output vector for the solution.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient, typename ASize>
	void solveInPlace(Function&& f, Gradient&& g, const ASize size, Scalar* x) requires (isSize<ASize>) { solveInPlace(makeOracle<Scalar>(size, std::forward<Function>(f), std::forward<Gradient>(g)), x); }

	/**
	 * @brief Solve in place using raw function and gradient functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param x Initial guess and output vector for the solution.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient>
	void solveInPlace(Function&& f, Gradient&& g, const std::span<Scalar> x) { solveInPlace(makeOracle<Scalar>(Size(x.size()), std::forward<Function>(f), std::forward<Gradient>(g)), x); }
	
	/**
	 * @brief Solve with initial guess using raw function and gradient functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param x0 Initial guess.
	 * @param size Problem dimension.
	 * @param x Output solution vector.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient, typename ASize>
	void solveWithGuess(Function&& f, Gradient&& g, const Scalar* x0, const ASize size, Scalar* x) requires (isSize<ASize>) { solveWithGuess(makeOracle<Scalar>(size, std::forward<Function>(f), std::forward<Gradient>(g)), x0, x); }

	/**
	 * @brief Solve with initial guess using raw function and gradient functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param x0 Initial guess.
	 * @param x Output solution vector.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient, std::ranges::input_range InitialGuess>
	void solveWithGuess(Function&& f, Gradient&& g, const InitialGuess&& x0, const std::span<Scalar> x) { solveWithGuess(makeOracle<Scalar>(Size(x.size()), std::forward<Function>(f), std::forward<Gradient>(g)), std::forward<InitialGuess>(x0), x); }
	
	// ========================================================================
	// MONITORING METHODS
	// ========================================================================
	
	Scalar getValue      () const { return m_fx;       } ///<  @brief Get last computed function value.
	Scalar getError      () const { return m_residual; } ///<  @brief Get last computed gradient norm.
	Size   getMaxIt      () const { return m_maxIt;    } ///<  @brief Get the maximum number of iterations allowed.
	Scalar getRelTol     () const { return m_relTol;   } ///<  @brief Get the convergence relative tolerance.
	Scalar getAbsTol     () const { return m_absTol;   } ///<  @brief Get the convergence absolute tolerance.
	Size   getIterations () const { return m_nIt;      } ///<  @brief Get the actual number of iterations performed.
	Info   getInfo       () const { return m_info;     } ///<  @brief Get the solver exit status.
	
	Size getInnerIterations(const Size it) const { return m_innerIts[it]; } ///<  @brief  Get number of inner iterations required to compute a step.
	
	const Criterion& getCriterion() const { return m_criterion; }
	
	Criterion& getCriterion() { return m_criterion; }
	
	void setMaxIt  (const Size    maxIt) { m_maxIt  = maxIt; } ///<  @brief Set the maximum number of iterations.
	void setRelTol (const Scalar& tol)   { m_relTol = tol;   } ///<  @brief Set the convergence relative tolerance.
	void setAbsTol (const Scalar& tol)   { m_absTol = tol;   } ///<  @brief Set the convergence absolute tolerance.
	
	void setTol(const Scalar& tol) { setRelTol(tol); setAbsTol(tol); } ///<  @brief Set both relative and absolute tolerance.
	
	/**
	 * @brief Set a file stream for solver logging output (e.g., stdout or file).
	 * @param out File pointer for output (nullptr disables output).
	 */
	void setOutput(std::FILE* out) { m_out = out; }
protected:
	Size   m_maxIt;    ///<  @brief Maximum number of iterations
	Scalar m_relTol;   ///<  @brief Relative tolerance for convergence
	Scalar m_absTol;   ///<  @brief Absolute tolerance for convergence
	Size   m_nIt;      ///<  @brief Number of iterations actually performed
	Info   m_info;     ///<  @brief Status of the solver after termination
	Scalar m_fx;       ///<  @brief Function value at last iteration
	Scalar m_residual; ///<  @brief A norm of the gradient at last iteration
	
	Criterion m_criterion; ///<  @brief Convergence criterion for the solver
	
	Size m_workCapacity = Size{}; ///<  @brief Maximum size of the arrays allocated by the solver
	
	std::vector< Size > m_innerIts; ///<  @brief Inner iterations tracking
	
	std::FILE* m_out = nullptr; ///<  @brief Optional output stream
};

#define LNOT_DEFINE_FIRST_ORDER_SOLVER \
	using Base   = FirstOrderSolverBase<Self>; \
	using Size   = typename Base::Size; \
	using Scalar = typename Base::Scalar; \
	using Info   = typename Base::Info; \
	\
	template<typename ASize> static constexpr bool isSize = Base::template isSize<ASize>; \
	template<typename ABool> static constexpr bool isBool = Base::template isBool<ABool>; \

#define LNOT_FIRST_ORDER_SOLVER_ATTRIBUTE \
	using Base::m_maxIt; \
	using Base::m_relTol; \
	using Base::m_absTol; \
	using Base::m_nIt; \
	using Base::m_info; \
	using Base::m_fx; \
	using Base::m_residual; \
	using Base::m_criterion; \
	using Base::m_workCapacity; \
	using Base::m_innerIts; \
	using Base::m_out; \

template<class Solver> concept CFirstOrderSolver = std::derived_from<Solver, FirstOrderSolverBase<Solver>>;

} // namespace LightOptim

#endif // LNOT_FIRST_ORDER_SOLVER_BASE_HPP
