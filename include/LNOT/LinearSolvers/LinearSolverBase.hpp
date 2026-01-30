#ifndef LNOT_LINEAR_SOLVER_BASE_HPP
#define LNOT_LINEAR_SOLVER_BASE_HPP

#include <cstdio>

#include <type_traits>
#include <limits>
#include <algorithm>

#include <LNOT/BasicLinalg/IdentityPreconditionerOp.hpp>
#include <LNOT/CRTPBase.hpp>

#include <BIC/Core.hpp>

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
 * - solveImpl() with an without a preconditioner
 */
template<class Derived>
class LinearSolverBase : public CRTPBase<Derived>
{
	using DerivedTraits = LinearSolverTraits<Derived>;
	using CRTP          = CRTPBase<Derived>;
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
	
	using IdOp = IdentityPreconditionerOp<Scalar>;
	
	template<typename HesOp>                  struct IsHessianOp   : BIC::Fixed<bool, std::is_invocable<HesOp, const Scalar*, Scalar*>::value > {}; ///<  @brief Trait to check if a type is a valid Hessian operator.
	template<typename HesOp, typename PrecOp> struct AreHessianOps : BIC::Fixed<bool, IsHessianOp<HesOp>::value and IsHessianOp<PrecOp>::value> {}; ///<  @brief Trait to check if two types are both valid Hessian operators.
	template<typename ASize>                  struct IsSize        : BIC::Fixed<bool, std::is_same<Size, BIC::Mutable<ASize>>::value > {};          ///<  @brief Trait to check if a type is either a `Size` or a `BIC::Fixed<Size, VALUE>`
	
	static constexpr Scalar defaultEps = std::numeric_limits<Scalar>::epsilon(); ///<  @brief Default value for relative and absolute tolerance of the solver. 
	
	/**
	 * @brief Construct the solver with optional maximum iterations and tolerance.
	 * @param maxIt Maximum number of iterations (default: 200000).
	 * @param relTol Convergence relative tolerance (default: machine epsilon).
	 * @param absTol Convergence absolute tolerance (default: machine epsilon).
	 */
	LinearSolverBase(const Size maxIt = 200000, const Scalar relTol = defaultEps, const Scalar absTol = defaultEps) : m_maxIt(maxIt), m_relTol(relTol), m_absTol(absTol) {}
	~LinearSolverBase() { clearWorkSpace(); }
	
	void clearWorkSpace() { CRTP::derived().clearWorkSpace(); } ///<  @brief Clear any internal memory or workspace used by the solver.
	
	/**
	 * @brief Solve the linear system Hx = -g using the provided Hessian operator.
	 * @param H Hessian operator (must satisfy IsHessianOp).
	 * @param g Gradient vector.
	 * @param size Problem size (must satisfy IsSize).
	 * @param x Solution vector (output).
	 */
	template<typename Op, typename ASize> 
	void solve(const Op& H, const Scalar* g, const ASize size, Scalar* x) requires (IsHessianOp<Op>::value and IsSize<ASize>::value) { IdOp I(size); solve(H, I, g, size, x); }
	
	/**
	 * @brief Solve the preconditioned linear system Hx = -g using the provided Hessian operator.
	 * @param H Hessian operator (must satisfy IsHessianOp).
	 * @param invB a preconditioner operator (must satisfy IsHessianOp).
	 * @param g Gradient vector.
	 * @param size Problem size (must satisfy IsSize).
	 * @param x Solution vector (output).
	 */
	template<typename HesOp, typename PrecOp, typename ASize> 
	void solve(const HesOp& H, const PrecOp& invB, const Scalar* g, const ASize size, Scalar* x) requires (AreHessianOps<HesOp,PrecOp>::value and IsSize<ASize>::value) { CRTP::derived().solveImpl(H, invB, g, size, BIC::fixed<bool, false>, x); }
	
	/**
	 * @brief Solve the system with an initial guess.
	 * @param H Hessian operator  (must satisfy IsHessianOp).
	 * @param g Gradient vector.
	 * @param x0 Initial guess vector.
	 * @param size Problem size (must satisfy IsSize).
	 * @param x Solution vector (output).
	 */
	template<typename Op, typename ASize> 
	void solveWithGuess(const Op& H, const Scalar* g, const Scalar* x0, const ASize size, Scalar* x) requires (IsHessianOp<Op>::value and IsSize<ASize>::value) { IdOp I(size); solveWithGuess(H, I, g, x0, size, x);  }
	
	/**
	 * @brief Solve the preconditioned system with an initial guess.
	 * @param H Hessian operator.
	 * @param invB a preconditioner operator (must satisfy IsHessianOp).
	 * @param g Gradient vector.
	 * @param x0 Initial guess vector.
	 * @param size Problem size (must satisfy IsSize).
	 * @param x Solution vector (output).
	 */	
	template<typename HesOp, typename PrecOp, typename ASize>  
	void solveWithGuess(const HesOp& H, const PrecOp& invB, const Scalar* g, const Scalar* x0, const ASize size, Scalar* x) requires (AreHessianOps<HesOp,PrecOp>::value and IsSize<ASize>::value) { std::copy(x0, x0 + size, x); CRTP::derived().solveImpl(H, invB, g, size, BIC::fixed<bool, true>, x);  }
	
	Scalar getError        () const { return CRTP::derived().getErrorImpl();        } ///<  @brief Get the final error after solving. Delegates to `Derived::getError()`.
	Scalar getSquaredError () const { return CRTP::derived().getSquaredErrorImpl(); } ///<  @brief Get the final squared error after solving. Delegates to `Derived::getSquaredError()`.
	
	Size   getMaxIt      () const { return m_maxIt;  } ///<  @brief Get the maximum number of iterations allowed
	Scalar getRelTol     () const { return m_relTol; } ///<  @brief Get the convergence relative tolerance.
	Scalar getAbsTol     () const { return m_absTol; } ///<  @brief Get the convergence absolute tolerance.
	Size   getIterations () const { return m_nIt;    } ///<  @brief Get the actual number of iterations performed.
	Info   getInfo       () const { return m_info;   } ///<  @brief Get the solver exit status.
	
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
	Size   m_maxIt;   ///<  @brief Maximum number of iterations
	Scalar m_relTol;  ///<  @brief Relative tolerance for convergence
	Scalar m_absTol;  ///<  @brief Absolute tolerance for convergence
	Size   m_nIt;     ///<  @brief Number of iterations actually performed
	Info   m_info;    ///<  @brief Status of the solver after termination
	
	Size m_workCapacity = 0; ///<  @brief Maximum size of the arrays allocated by the solver
	
	std::FILE* m_out = nullptr; ///<  @brief Optional output stream
};

#define LNOT_DEFINE_LINEAR_SOLVER \
	using Base   = LinearSolverBase<Self>; \
	using Size   = typename Base::Size; \
	using Scalar = typename Base::Scalar; \
	using Info   = typename Base::Info; \
	\
	template<typename HesOp, typename PrecOp> using AreHessianOps = typename Base::template AreHessianOps<HesOp,PrecOp>; \
	template<typename ASize>                  using IsSize        = typename Base::template IsSize<ASize>; \
	
#define LNOT_LINEAR_SOLVER_ATTRIBUTE \
	using Base::m_maxIt; \
	using Base::m_relTol; \
	using Base::m_absTol; \
	using Base::m_nIt; \
	using Base::m_info; \
	using Base::m_workCapacity; \
	using Base::m_out; \

template<class T> struct IsLinearSolver : BIC::Fixed<bool, std::is_base_of<LinearSolverBase<T>, T>::value > {};

template<class T> concept CLinearSolver = IsLinearSolver<T>::value;

} // namespace LNOT

#endif // LNOT_LINEAR_SOLVER_BASE_HPP
