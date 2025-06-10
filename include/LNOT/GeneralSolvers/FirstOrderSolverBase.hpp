#ifndef LNOT_FIRST_ORDER_MINIMIZER_BASE_HPP
#define LNOT_FIRST_ORDER_MINIMIZER_BASE_HPP

#include <LNOT/Oracles/OracleWrapper.hpp>

#include <cstdio>
#include <cmath>

#include <vector>
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
class FirstOrderSolverBase
{
	using DerivedTraits = FirstOrderSolverTraits<Derived>;
public:
	using Scalar = typename DerivedTraits::Scalar; ///<  @brief The scalar type used in computations (e.g., float, double)
	using Size   = typename DerivedTraits::Size;   ///<  @brief The size type used for indexing and loop counters
	
	/// @brief Enumeration indicating solver termination status.
	enum class Info 
	{
		SUCCESS,  ///<  Converged successfully
		FAILURE,  ///<  Failed to converge
		BREAKDOWN ///<  Numerical breakdown (e.g., division by zero)
	};
	
	const Derived& derived_cast() const { return static_cast<const Derived&>(*this); }
	      Derived& derived_cast()       { return static_cast<      Derived&>(*this); }
	
	template<typename Function>  struct IsFunction  : std::bool_constant< std::is_invocable<Function, const Scalar*>::value > {};                           ///<  @brief Trait to check if Function is a valid callable with signature Scalar(const Scalar*)
	template<typename Gradient>  struct IsGradient  : std::bool_constant< std::is_invocable<Gradient, const Scalar*, Scalar*>::value > {};                  ///<  @brief Trait to check if Gradient is a valid callable with signature void(const Scalar*, Scalar*)
	template<typename Function, typename Gradient> struct IsProgram : std::bool_constant< IsFunction<Function>::value and IsGradient<Gradient>::value > {}; ///<  @brief Trait to check if Function + Gradient pair is valid
	
	// ========================================================================
	// SOLVER INTERFACES
	// ========================================================================
	
	/**
	 * @brief Constructor
	 * @param maxIt Maximum number of iterations (default: 200000).
	 * @param tol Convergence tolerance (default: machine epsilon).
	 */
	FirstOrderSolverBase(const Size maxIt = 200000, const Scalar tol = std::numeric_limits<Scalar>::epsilon()) : m_maxIt(maxIt), m_tol(tol) {}
	~FirstOrderSolverBase() { clearWorkSpace(); }
	
	/// @brief Clear any internal memory or workspace used by the solver.
	void clearWorkSpace() { derived_cast().clearWorkSpace(); }
		
	/**
	 * @brief Solve using a valid FirstOrderOracle (no initial guess).
	 * @param oracle An oracle object conforming to FirstOrderOracle_concept.
	 * @param x Output vector for the solution.
	 */
	template<FirstOrderOracle_concept Oracle>
	void solve(Oracle& oracle, Scalar* x) { solve_impl(oracle, std::false_type{}, x); }
	
	/**
	 * @brief Solve using a valid FirstOrderOracle with an initial guess.
	 * @param oracle An oracle object.
	 * @param x0 Initial guess.
	 * @param x Output vector for the solution.
	 */
	template<FirstOrderOracle_concept Oracle>
	void solveWithGuess(Oracle& oracle, const Scalar* x0, Scalar* x) { std::copy(x0, x0 + oracle.getNDims(), x); solve_impl(oracle, std::true_type{}, x); }
	
	/**
	 * @brief Solve using raw function and gradient functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param size Problem dimension.
	 * @param x Output solution vector.
	 */
	template<typename Function, typename Gradient>
	void solve(Function f, Gradient g, const Size size, Scalar* x) requires (IsProgram<Function,Gradient>::value) { solve_impl(f, g, size, std::false_type{}, x); }
	
	/**
	 * @brief Solve with initial guess using raw function and gradient functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param x0 Initial guess.
	 * @param size Problem dimension.
	 * @param x Output solution vector.
	 */
	template<typename Function, typename Gradient>
	void solveWithGuess(Function f, Gradient g, const Scalar* x0, const Size size, Scalar* x) requires (IsProgram<Function,Gradient>::value) { std::copy(x0, x0 + size, x); solve_impl(f, g, size, std::true_type{}, x);  }
	
	/**
	 * @brief Internal function to create an OracleWrapper from functors and solve with or without an initial guess.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param size Problem dimension.
	 * @param bc a placeholder for `solveInPlace`.
	 * @param x Output solution vector.
	 * 
	 * @tparam solveInPlace specifying if x should be used as an initial guess.
	 */
	template<typename Function, typename Gradient, bool solveInPlace> 
	void solve_impl(Function f, Gradient g, const Size size, std::bool_constant<solveInPlace> bc, Scalar* x) requires (IsProgram<Function,Gradient>::value) { OracleWrapper<Scalar,Function,Gradient> oracle(size, f, g); solve_impl(oracle, bc, x); }
	
	/**
	 * @brief Solve using a valid FirstOrderOracle with or without an initial guess.
	 * @param oracle An oracle object.
	 * @param bc a placeholder for `solveInPlace`.
	 * @param x Output solution vector.
	 * 
	 * Delegates to `Derived::solve_impl()`.
	 * 
	 * @tparam solveInPlace specifying if x should be used as an initial guess.
	 */
	template<FirstOrderOracle_concept Oracle, bool solveInPlace> 
	void solve_impl(Oracle& oracle, std::bool_constant<solveInPlace> bc, Scalar* x) { derived_cast().solve_impl(oracle, bc, x); }
	
	// ========================================================================
	// MONITORING METHODS
	// ========================================================================
	
	Scalar getValue        () const { return m_fx; }                         ///<  @brief Get last computed function value.
	Scalar getError        () const { return std::sqrt(m_squaredNormGrad); } ///<  @brief Get L2 norm of the last computed gradient.
	Scalar getSquaredError () const { return m_squaredNormGrad; }            ///<  @brief Get squared L2 norm of the last computed gradient.
	
	Size   getMaxIt      () const { return m_maxIt; } ///<  @brief Get the maximum number of iterations allowed.
	Scalar getTol        () const { return m_tol;   } ///<  @brief Get the convergence tolerance.
	Size   getIterations () const { return m_nIt;   } ///<  @brief Get the actual number of iterations performed.
	Info   getInfo       () const { return m_info;  } ///<  @brief Get the solver exit status.
	
	Size getInnerIterations(const Size it) const { return m_innerIts[it]; } ///<  @brief  Get number of inner iterations required to compute a step.
	
	void setMaxIt (const Size    maxIt) { m_maxIt = maxIt; } ///<  @brief Set the maximum number of iterations.
	void setTol   (const Scalar& tol)   { m_tol   = tol;   } ///<  @brief Set the convergence tolerance.
	
	/**
	 * @brief Set a file stream for solver logging output (e.g., stdout or file).
	 * @param out File pointer for output (nullptr disables output).
	 */
	void setOutput(std::FILE* out) { m_out = out; }
protected:
	Size   m_maxIt;           ///<  @brief Maximum number of iterations
	Scalar m_tol;             ///<  @brief Tolerance for convergence
	Size   m_nIt;             ///<  @brief Number of iterations actually performed
	Info   m_info;            ///<  @brief Status of the solver after termination
	Scalar m_squaredNormGrad; ///<  @brief Squared norm of the gradient
	Scalar m_fx;              ///<  @brief Function value at last iteration
	
	Size m_workCapacity = 0; ///<  @brief Maximum size of the arrays allocated by the solver
	
	std::vector< Size > m_innerIts; ///<  @brief Inner iterations tracking
	
	std::FILE* m_out = nullptr; ///<  @brief Optional output stream
};

template<class T> struct IsFirstOrderSolver : std::bool_constant< std::is_base_of<FirstOrderSolverBase<T>, T>::value > {}; ///<  @brief Trait to determine if a type derives from FirstOrderSolverBase.

template<class T> concept IsFirstOrderSolver_concept = IsFirstOrderSolver<T>::value;

} // namespace LightOptim

#endif // LNOT_FIRST_ORDER_MINIMIZER_BASE_HPP
