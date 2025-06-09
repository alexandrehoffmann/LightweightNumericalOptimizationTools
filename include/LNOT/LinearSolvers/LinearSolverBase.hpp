#ifndef LNOT_LINEAR_SOLVER_BASE_HPP
#define LNOT_LINEAR_SOLVER_BASE_HPP

#include <cstdio>

#include <type_traits>
#include <limits>
#include <algorithm>

namespace LNOT
{

template<class LinearSolver> struct LinearSolverTraits;

/**
 * @brief Base class for linear solvers using the CRTP (Curiously Recurring Template Pattern).
 * 
 * This class provides a generic interface for iterative linear solvers, where the actual implementation
 * is expected to be provided by the derived class. It supports both simple Hessian-vector operations and
 * preconditioned variants. It also tracks solver metadata such as iteration count and error.
 * 
 * @tparam Derived The derived class implementing specific solver behavior.
 * Derived must implement:
 * - void clearWorkSpace()
 * - solve_impl() with an without a preconditioner
 */
template<class Derived>
class LinearSolverBase
{
	using DerivedTraits = LinearSolverTraits<Derived>;
public:
	using Scalar = typename DerivedTraits::Scalar; ///<  @brief The scalar type used in computations (e.g., float, double)
	using Size   = typename DerivedTraits::Size;   ///<  @brief The size type used for indexing and loop counters
	/// @brief Enumeration indicating solver termination status.
	enum class Info 
	{
		SUCCESS,            ///<  Converged successfully
		FAILURE,            ///<  Failed to converge
		NEGATIVE_CURVATURE, ///<  Encountered negative curvature (The problem is ill-paused)
		BREAKDOWN           ///<  Numerical breakdown (e.g., division by zero)
	};
	
	template<typename HesOp>                  struct IsHessianOp   : std::bool_constant< std::is_invocable<HesOp, const Scalar*, Scalar*>::value > {}; ///<  @brief Trait to check if a type is a valid Hessian operator.
	template<typename HesOp, typename PrecOp> struct AreHessianOps : std::bool_constant< IsHessianOp<HesOp>::value and IsHessianOp<PrecOp>::value> {}; ///<  @brief Trait to check if two types are both valid Hessian operators.

	const Derived& derived_cast() const { return static_cast<const Derived&>(*this); }
	      Derived& derived_cast()       { return static_cast<      Derived&>(*this); }
	/**
	 * @brief Construct the solver with optional maximum iterations and tolerance.
	 * @param maxIt Maximum number of iterations (default: 200000).
	 * @param tol Convergence tolerance (default: machine epsilon).
	 */
	LinearSolverBase(const Size maxIt = 200000, const Scalar tol = std::numeric_limits<Scalar>::epsilon()) : m_maxIt(maxIt), m_tol(tol) {}
	~LinearSolverBase() { clearWorkSpace(); }
	/// @brief Clear any internal memory or workspace used by the solver.
	void clearWorkSpace() { derived_cast().clearWorkSpace(); }
	/**
	 * @brief Solve the linear system Hx = -g using the provided Hessian operator.
	 * @param H Hessian operator (must satisfy IsHessianOp).
	 * @param g Gradient vector.
	 * @param size Problem size (dimension).
	 * @param x Solution vector (output).
	 */
	template<typename Op> 
	void solve(const Op& H, const Scalar* g, const Size size, Scalar* x) requires (IsHessianOp<Op>::value) { solve_impl(H, g, size, std::false_type{}, x); }
	/**
	 * @brief Solve the system with an initial guess.
	 * @param H Hessian operator  (must satisfy IsHessianOp).
	 * @param g Gradient vector.
	 * @param x0 Initial guess vector.
	 * @param size Problem size.
	 * @param x Solution vector (output).
	 */
	template<typename Op> 
	void solveWithGuess(const Op& H, const Scalar* g, const Scalar* x0, const Size size, Scalar* x) requires (IsHessianOp<Op>::value) { std::copy(x0, x0 + size, x); solve_impl(H, g, size, std::true_type{}, x);  }
	/**
	 * @brief Solve the system with or without an initial guess.
	 * 
	 * Delegates to `Derived::solve_impl()`.
	 * 
	 * @param H Hessian operator (must satisfy IsHessianOp).
	 * @param g Gradient vector.
	 * @param size Problem size.
	 * @param x Solution vector (output if solveInPlace == false, input and output else).
	 * 
	 * @tparam solveInPlace specifying if x should be used as an initial guess.
	 */
	template<typename Op, bool solveInPlace> 
	void solve_impl(const Op& H, const Scalar* g, const Size size, std::bool_constant<solveInPlace> bc, Scalar* x) requires (IsHessianOp<Op>::value) { derived_cast().solve_impl(H, g, size, bc, x); }
	/**
	 * @brief Solve the preconditioned linear system Hx = -g using the provided Hessian operator.
	 * @param H Hessian operator (must satisfy IsHessianOp).
	 * @param invB a preconditioner operator (must satisfy IsHessianOp).
	 * @param g Gradient vector.
	 * @param size Problem size (dimension).
	 * @param x Solution vector (output).
	 */
	template<typename HesOp, typename PrecOp> 
	void solve(const HesOp& H, const PrecOp& invB, const Scalar* g, const Size size, Scalar* x) requires (AreHessianOps<HesOp,PrecOp>::value) { solve_impl(H, invB, g, size, std::false_type{}, x); }
	/**
	 * @brief Solve the preconditioned system with an initial guess.
	 * @param H Hessian operator.
	 * @param invB a preconditioner operator (must satisfy IsHessianOp).
	 * @param g Gradient vector.
	 * @param x0 Initial guess vector.
	 * @param size Problem size.
	 * @param x Solution vector (output).
	 */	
	template<typename HesOp, typename PrecOp> 
	void solveWithGuess(const HesOp& H, const PrecOp& invB, const Scalar* g, const Scalar* x0, const Size size, Scalar* x) requires (AreHessianOps<HesOp,PrecOp>::value) { std::copy(x0, x0 + size, x); solve_impl(H, invB, g, size, std::true_type{}, x);  }
	/**
	 * @brief Solve the preconditioned system with or without an initial guess.
	 * 
	 * Delegates to `Derived::solve_impl()`.
	 * 
	 * @param H Hessian operator (must satisfy IsHessianOp).
	 * @param invB a preconditioner operator (must satisfy IsHessianOp).
	 * @param g Gradient vector.
	 * @param size Problem size.
	 * @param x Solution vector (output if solveInPlace == false, input and output else).
	 * 
	 * @tparam solveInPlace specifying if x should be used as an initial guess.
	 */
	template<typename HesOp, typename PrecOp, bool solveInPlace> 
	void solve_impl(const HesOp& H, const PrecOp& invB, const Scalar* g, const Size size, std::bool_constant<solveInPlace> bc, Scalar* x) requires (AreHessianOps<HesOp,PrecOp>::value) { derived_cast().solve_impl(H, invB, g, size, bc, x); }
	
	Scalar getError        () const { return derived_cast().getError();        } ///<  @brief Get the final error after solving. Delegates to `Derived::getError()`.
	Scalar getSquaredError () const { return derived_cast().getSquaredError(); } ///<  @brief Get the final squared error after solving. Delegates to `Derived::getSquaredError()`.
	
	Size   getMaxIt      () const { return m_maxIt; } ///<  @brief Get the maximum number of iterations allowed
	Scalar getTol        () const { return m_tol;   } ///<  @brief Get the convergence tolerance.
	Size   getIterations () const { return m_nIt;   } ///<  @brief Get the actual number of iterations performed.
	Info   getInfo       () const { return m_info;  } ///<  @brief Get the solver exit status.
	
	void setMaxIt (const Size    maxIt) { m_maxIt = maxIt; } ///<  @brief Set the maximum number of iterations.
	void setTol   (const Scalar& tol)   { m_tol   = tol;   } ///<  @brief Set the convergence tolerance.
	
	/**
	 * @brief Set a file stream for solver logging output (e.g., stdout or file).
	 * @param out File pointer for output (nullptr disables output).
	 */
	void setOutput(std::FILE* out) { m_out = out; }
protected:
	Scalar getResidualThreshold()        const { return m_tol*std::max(Scalar(1), getError()); }              ///<  @brief Get the residual threshold for which the problem is considered solved.
	Scalar getSquaredResidualThreshold() const { return m_tol*m_tol*std::max(Scalar(1), getSquaredError()); } ///<  @brief Get the squared residual threshold for which the problem is considered solved. 

	Size   m_maxIt; ///<  Maximum number of iterations
	Scalar m_tol;   ///<  Tolerance for convergence
	Size   m_nIt;   ///<  Number of iterations actually performed
	Info   m_info;  ///<  Status of the solver after termination
	
	Size m_workCapacity = 0; ///<  Maximum size of the arrays allocated by the solver
	
	std::FILE* m_out = nullptr; ///<  Optional output stream
};

template<class T> struct IsLinearSolver : std::bool_constant< std::is_base_of<LinearSolverBase<T>, T>::value > {};

template<class T> concept LinearSolver_concept = IsLinearSolver<T>::value;

} // namespace LNOT

#endif // LNOT_LINEAR_SOLVER_BASE_HPP
