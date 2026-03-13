#ifndef LNOT_SECOND_ORDER_SOLVER_BASE_HPP
#define LNOT_SECOND_ORDER_SOLVER_BASE_HPP

#include <LNOT/Oracles/OracleWrapper.hpp>
#include <LNOT/Oracles/concepts.hpp>
#include <LNOT/FloatingPoint/NumTraits.hpp>
#include <BIC/Core.hpp>

#include <cstdio>
#include <cmath>

#include <vector>
#include <algorithm>

namespace LNOT
{

template<class Derived> struct SecondOrderSolverTraits;

/**
 * @brief Base class for second-order optimization solvers using the CRTP pattern.
 * 
 * This class provides generic solve interfaces for different kinds of user input:
 * raw function/gradient/Hessian functors, or a conforming oracle. It handles workspace
 * management, iteration control, output, and error tracking.
 *
 * @tparam Derived The derived solver class implementing the actual algorithm.
 */
template<class Derived>
class SecondOrderSolverBase : public CRTPBase<Derived>
{
	using DerivedTraits = SecondOrderSolverTraits<Derived>;
	using CRTP          = CRTPBase<Derived>;
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
	
	template<typename ABool> struct IsBool : BIC::Fixed<bool, std::is_same<bool, BIC::Mutable<ABool> >::value > {}; ///<  @brief Trait to check if a type is either a `bool` or a `BIC::Fixed<bool, VALUE>`	
	template<typename ASize> struct IsSize : BIC::Fixed<bool, std::is_same<Size, BIC::Mutable<ASize>>::value > {};  ///<  @brief Trait to check if a type is either a `Size` or a `BIC::Fixed<Size, VALUE>`	
	
	static inline const Scalar defaultEps = NumTraits<Scalar>::epsilon; ///<  @brief Default value for relative and absolute tolerance of the solver.  
	
	// ===================================================================
	// SOLVER INTERFACES
	// ===================================================================
	
	/**
	 * @brief Constructor
	 * @param maxIt Maximum number of iterations (default: 200000).
	 * @param tol Convergence tolerance (default: machine epsilon).
	 */
	SecondOrderSolverBase(const Size maxIt = 200000, const Scalar relTol = defaultEps, const Scalar absTol = defaultEps) : m_maxIt(maxIt), m_relTol(relTol), m_absTol(absTol) {}
	~SecondOrderSolverBase() { clearWorkSpace(); }
	
	/// @brief Clear any internal memory or workspace used by the solver.
	void clearWorkSpace() { CRTP::derived().clearWorkSpaceImpl(); }
	
	/**
	 * @brief Solve using a valid SecondOrderOracle (no initial guess).
	 * @param oracle An oracle object conforming to CSecondOrderOracle.
	 * @param x Output vector for the solution.
	 */
	template<CSecondOrderOracle Oracle>
	void solve(Oracle& oracle, Scalar* x) { CRTP::derived().solveImpl(oracle, BIC::fixed<bool,true>, x); }
	
	/**
	 * @brief Solve using a valid SecondOrderOracle with an initial guess.
	 * @param oracle An oracle object.
	 * @param x0 Initial guess.
	 * @param x Output vector for the solution.
	 */
	template<CSecondOrderOracle Oracle>
	void solveWithGuess(Oracle& oracle, const Scalar* x0, Scalar* x) { std::copy(x0, x0 + oracle.getNDims(), x); CRTP::derived().solveImpl(oracle, BIC::fixed<bool,false>, x); }
	
	/**
	 * @brief Solve using raw function, gradient and Hessian product functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param H Hessian-product functor.
	 * @param size Problem dimension.
	 * @param x Output solution vector.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient, CHessianOp<Scalar> HessianOp, typename ASize>
	void solve(Function f, Gradient g, HessianOp H, const ASize size, Scalar* x) requires (IsSize<ASize>::value) { solve(f, g, H, size, BIC::fixed<bool,false>, x); }
	
	/**
	 * @brief Solve with initial guess using raw function, gradient and Hessian product functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param H Hessian-product functor.
	 * @param x0 Initial guess.
	 * @param size Problem dimension.
	 * @param x Output solution vector.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient, CHessianOp<Scalar> HessianOp, typename ASize>
	void solveWithGuess(Function f, Gradient g, HessianOp H, const Scalar* x0, const ASize size, Scalar* x) requires (IsSize<ASize>::value) { std::copy(x0, x0 + size, x); CRTP::derived().solveImpl(f, g, H, size, BIC::fixed<bool,true>, x);  }
	
	/**
	 * @brief Solve using raw function, gradient, Hessian product and preconditioner operator functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param H Hessian-product functor.
	 * @param invB Preconditioner functor.
	 * @param size Problem dimension.
	 * @param x Output solution vector.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient, CHessianOp<Scalar> HessianOp, CHessianOp<Scalar> PrecondOp, typename ASize>
	void solve(Function f, Gradient g, HessianOp H, PrecondOp invB, const ASize size, Scalar* x) requires (IsSize<ASize>::value) { solve(f, g, H, invB, size, BIC::fixed<bool,false>, x); }
	
	/**
	 * @brief Solve with an initial guess using raw function, gradient, Hessian product and preconditioner operator functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param H Hessian-product functor.
	 * @param invB Preconditioner functor.
	 * @param x0 Initial guess.
	 * @param size Problem dimension.
	 * @param x Output solution vector.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient, CHessianOp<Scalar> HessianOp, CHessianOp<Scalar> PrecondOp, typename ASize>
	void solveWithGuess(Function f, Gradient g, HessianOp H, PrecondOp invB, const Scalar* x0, const ASize size, Scalar* x) requires (IsSize<ASize>::value) { std::copy(x0, x0 + size, x); solve(f, g, H, invB, size, BIC::fixed<bool,true>, x);  }
	
	/**
	 * @brief Internal function to create an OracleWrapper from functors and solve with or without an initial guess.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param H Hessian-product functor.
	 * @param size Problem dimension.
	 * @param bc a placeholder for `solveInPlace`.
	 * @param x Output solution vector.
	 * 
	 * @tparam solveInPlace specifying if x should be used as an initial guess.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient, CHessianOp<Scalar> HessianOp, typename ASize, typename ABool> 
	void solve(Function f, Gradient g, HessianOp H, const ASize size, const ABool solveInPlace, Scalar* x) requires (IsSize<ASize>::value and IsBool<ABool>::value) { OracleWrapper<Scalar,ASize,Function,Gradient,HessianOp> oracle(size, f, g, H); CRTP::derived().solveImpl(oracle, solveInPlace, x); }
	
	/**
	 * @brief Internal function to create an OracleWrapper from functors and solve with or without an initial guess.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param H Hessian-product functor.
	 * @param invB Preconditioner functor.
	 * @param size Problem dimension.
	 * @param bc a placeholder for `solveInPlace`.
	 * @param x Output solution vector.
	 * 
	 * @tparam solveInPlace specifying if x should be used as an initial guess.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient, CHessianOp<Scalar> HessianOp, CHessianOp<Scalar> PrecondOp, typename ASize, typename ABool> 
	void solve(Function f, Gradient g, HessianOp H, PrecondOp invB, const ASize size, const ABool solveInPlace, Scalar* x) requires (IsSize<ASize>::value and IsBool<ABool>::value) { OracleWrapper<Scalar,ASize,Function,Gradient,HessianOp,PrecondOp> oracle(size, f, g, H, invB); CRTP::derived().solveImpl(oracle, solveInPlace, x); }
	
	/**
	 * @brief Solve using a valid FirstOrderOracle with or without an initial guess.
	 * @param oracle An oracle object.
	 * @param bc a placeholder for `solveInPlace`.
	 * @param x Output solution vector.
	 * 
	 * Delegates to `Derived::solveImpl()`.
	 * 
	 * @tparam solveInPlace specifying if x should be used as an initial guess.
	 */
	template<CSecondOrderOracle Oracle, typename ABool> 
	void solve(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value) { CRTP::derived().solveImpl(oracle, solveInPlace, x); }
	
	// ===================================================================
	// MONITORING METHODS
	// ===================================================================
	
	Scalar getValue        () const { return m_fx; }                         ///<  @brief Get last computed function value
	Scalar getError        () const { return sqrt(m_squaredNormGrad); } ///<  @brief Get L2 norm of the last computed gradient
	Scalar getSquaredError () const { return m_squaredNormGrad; }            ///<  @brief Get squared L2 norm of the last computed gradient
	
	Size   getMaxIt      () const { return m_maxIt;  } ///<  @brief Get the maximum number of iterations allowed.
	Scalar getRelTol     () const { return m_relTol; } ///<  @brief Get the convergence relative tolerance.
	Scalar getAbsTol     () const { return m_absTol; } ///<  @brief Get the convergence absolute tolerance.
	Size   getIterations () const { return m_nIt;    } ///<  @brief Get the actual number of iterations performed.
	Info   getInfo       () const { return m_info;   } ///<  @brief Get the solver exit status.
	
	Size getInnerIterations(const Size it) const { return m_innerIts[it]; } ///<  @brief  Get number of inner iterations required to compute a step.
	
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
	Size   m_maxIt;           ///<  @brief Maximum number of iterations
	Scalar m_relTol;          ///<  @brief Relative tolerance for convergence
	Scalar m_absTol;          ///<  @brief Absolute tolerance for convergence
	Size   m_nIt;             ///<  @brief Number of iterations actually performed
	Info   m_info;            ///<  @brief Status of the solver after termination
	Scalar m_squaredNormGrad; ///<  @brief Squared norm of the gradient
	Scalar m_fx;              ///<  @brief Function value at last iteration
	
	Size m_workCapacity = Size{}; ///<  @brief Maximum size of the arrays allocated by the solver
	
	std::vector< Size > m_innerIts; ///<  @brief Inner iterations tracking
	
	std::FILE* m_out = nullptr; ///<  @brief Optional output stream
};

#define LNOT_DEFINE_SECOND_ORDER_SOLVER \
	using Base   = SecondOrderSolverBase<Self>; \
	using Size   = typename Base::Size; \
	using Scalar = typename Base::Scalar; \
	using Info   = typename Base::Info; \
	\
	template<typename ABool> using IsBool = typename Base::template IsBool<ABool>; \
	template<typename ASize> using IsSize = typename Base::template IsSize<ASize>; \

#define LNOT_SECOND_ORDER_SOLVER_ATTRIBUTE \
	using Base::m_maxIt; \
	using Base::m_relTol; \
	using Base::m_absTol; \
	using Base::m_nIt; \
	using Base::m_info; \
	using Base::m_squaredNormGrad; \
	using Base::m_fx; \
	using Base::m_workCapacity; \
	using Base::m_innerIts; \
	using Base::m_out; \

template<class T> struct IsSecondOrderSolver : std::bool_constant< std::is_base_of<SecondOrderSolverBase<T>, T>::value > {};  ///<  @brief Trait to determine if a type derives from FirstOrderSolverBase.

template<class T> concept CSecondOrderSolver = IsSecondOrderSolver<T>::value;

} // namespace LightOptim

#endif // LNOT_SECOND_ORDER_SOLVER_BASE_HPP
