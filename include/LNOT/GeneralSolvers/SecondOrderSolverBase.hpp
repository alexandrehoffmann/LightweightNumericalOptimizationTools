#ifndef LNOT_SECOND_ORDER_SOLVER_BASE_HPP
#define LNOT_SECOND_ORDER_SOLVER_BASE_HPP

#include <LNOT/Oracles/OracleWrapper.hpp>
#include <LNOT/Oracles/concepts.hpp>
#include <LNOT/FloatingPoint/NumTraits.hpp>
#include <BIC/Core.hpp>

#include <cstdio>
#include <cmath>

#include <vector>
#include <span>
#include <ranges>
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
	
	static inline constexpr Size order = 2;
	
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
	
	/**
	 * @brief Solve using a valid SecondOrderOracle (no initial guess).
	 * @param oracle An oracle object conforming to CSecondOrderOracle.
	 * @param x Output vector for the solution.
	 */
	template<CSecondOrderOracle Oracle>
	void solve(Oracle&& oracle, Scalar* x) { CRTP::derived().solveImpl(std::forward<Oracle>(oracle), BIC::fixed<bool,false>, x); }

	/**
	 * @brief Solve using a valid SecondOrderOracle (no initial guess).
	 * @param oracle An oracle object conforming to CSecondOrderOracle.
	 * @param x Output vector for the solution.
	 */
	template<CSecondOrderOracle Oracle>
	void solve(Oracle&& oracle, const std::span<Scalar> x) { assert(x.size() == oracle.getNDims()); solve(std::forward<Oracle>(oracle), x.data()); }
		
	/**
	 * @brief Solve using a valid SecondOrderOracle in place (uses x as an initial guess).
	 * @param oracle An oracle object conforming to CSecondOrderOracle.
	 * @param x Initial guess and output vector for the solution.
	 */
	template<CSecondOrderOracle Oracle>
	void solveInPlace(Oracle&& oracle, Scalar* x) { CRTP::derived().solveImpl(std::forward<Oracle>(oracle), BIC::fixed<bool,true>, x); }

	/**
	 * @brief Solve using a valid SecondOrderOracle in place (uses x as an initial guess).
	 * @param oracle An oracle object conforming to CSecondOrderOracle.
	 * @param x Initial guess and output vector for the solution.
	 */
	template<CSecondOrderOracle Oracle>
	void solveInPlace(Oracle&& oracle, const std::span<Scalar> x) { assert(x.size() == oracle.getNDims()); solveInPlace(std::forward<Oracle>(oracle), x.data()); }
	
	/**
	 * @brief Solve using a valid SecondOrderOracle with an initial guess.
	 * @param oracle An oracle object.
	 * @param x0 Initial guess.
	 * @param x Output vector for the solution.
	 */
	template<CSecondOrderOracle Oracle>
	void solveWithGuess(Oracle&& oracle, const Scalar* x0, Scalar* x) { std::copy(x0, x0 + oracle.getNDims(), x); solveInPlace(std::forward<Oracle>(oracle), x); }

	/**
	 * @brief Solve using a valid SecondOrderOracle with an initial guess.
	 * @param oracle An oracle object.
	 * @param x0 Initial guess.
	 * @param x Output vector for the solution.
	 */
	template<CSecondOrderOracle Oracle, std::ranges::input_range InitialGuess>
	void solveWithGuess(Oracle&& oracle, InitialGuess&& x0, const std::span<Scalar> x) { std::ranges::copy(x0, x.begin()); solveInPlace(std::forward<Oracle>(oracle), x); }
	
	/**
	 * @brief Solve using raw function, gradient and Hessian product functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param H Hessian-product functor.
	 * @param size Problem dimension.
	 * @param x Output solution vector.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient, CHessianOp<Scalar> HessianOp, typename ASize>
	void solve(Function&& f, Gradient&& g, HessianOp&& H, const ASize size, Scalar* x) requires (isSize<ASize>) { return solve(makeOracle<Scalar>(size, std::forward<Function>(f), std::forward<Gradient>(g), std::forward<HessianOp>(H)), x); }

	/**
	 * @brief Solve using raw function, gradient and Hessian product functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param H Hessian-product functor.
	 * @param x Output solution vector.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient, CHessianOp<Scalar> HessianOp>
	void solve(Function&& f, Gradient&& g, HessianOp&& H, const std::span<Scalar> x) { return solve(makeOracle<Scalar>(Size(x.size()), std::forward<Function>(f), std::forward<Gradient>(g), std::forward<HessianOp>(H)), x); }
	
	/**
	 * @brief Solve in place using raw function, gradient and Hessian product functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param H Hessian-product functor.
	 * @param size Problem dimension.
	 * @param x Initial guess and output vector for the solution.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient, CHessianOp<Scalar> HessianOp, typename ASize>
	void solveInPlace(Function&& f, Gradient&& g, HessianOp&& H, const ASize size, Scalar* x) requires (isSize<ASize>) { return solveInPlace(makeOracle<Scalar>(size, std::forward<Function>(f), std::forward<Gradient>(g), std::forward<HessianOp>(H)), x); }

	/**
	 * @brief Solve in place using raw function, gradient and Hessian product functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param H Hessian-product functor.
	 * @param x Initial guess and output vector for the solution.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient, CHessianOp<Scalar> HessianOp>
	void solveInPlace(Function&& f, Gradient&& g, HessianOp&& H, const std::span<Scalar> x) { return solveInPlace(makeOracle<Scalar>(Size(x.size()), std::forward<Function>(f), std::forward<Gradient>(g), std::forward<HessianOp>(H)), x); }
	
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
	void solveWithGuess(Function&& f, Gradient&& g, HessianOp&& H, const Scalar* x0, const ASize size, Scalar* x) requires (isSize<ASize>) { return solveWithGuess(makeOracle<Scalar>(size, std::forward<Function>(f), std::forward<Gradient>(g), std::forward<HessianOp>(H)), x0, x); }

	/**
	 * @brief Solve with initial guess using raw function, gradient and Hessian product functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param H Hessian-product functor.
	 * @param x0 Initial guess.
	 * @param x Output solution vector.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient, CHessianOp<Scalar> HessianOp, std::ranges::input_range InitialGuess>
	void solveWithGuess(Function&& f, Gradient&& g, HessianOp&& H, InitialGuess&& x0, const std::span<Scalar> x) { return solveWithGuess(makeOracle<Scalar>(Size(x.size()), std::forward<Function>(f), std::forward<Gradient>(g), std::forward<HessianOp>(H)), std::forward<InitialGuess>(x0), x); }
	
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
	void solve(Function&& f, Gradient&& g, HessianOp&& H, PrecondOp&& invB, const ASize size, Scalar* x) requires (isSize<ASize>) { return solve(makeOracle<Scalar>(size, std::forward<Function>(f), std::forward<Gradient>(g), std::forward<HessianOp>(H), std::forward<PrecondOp>(invB)), x); }

	/**
	 * @brief Solve using raw function, gradient, Hessian product and preconditioner operator functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param H Hessian-product functor.
	 * @param invB Preconditioner functor.
	 * @param x Output solution vector.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient, CHessianOp<Scalar> HessianOp, CHessianOp<Scalar> PrecondOp>
	void solve(Function&& f, Gradient&& g, HessianOp&& H, PrecondOp&& invB, const std::span<Scalar> x) { return solve(makeOracle<Scalar>(Size(x.size()), std::forward<Function>(f), std::forward<Gradient>(g), std::forward<HessianOp>(H), std::forward<PrecondOp>(invB)), x); }
	
	/**
	 * @brief Solve in place using raw function, gradient, Hessian product and preconditioner operator functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param H Hessian-product functor.
	 * @param invB Preconditioner functor.
	 * @param size Problem dimension.
	 * @param x Initial guess and output vector for the solution.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient, CHessianOp<Scalar> HessianOp, CHessianOp<Scalar> PrecondOp, typename ASize>
	void solveInPlace(Function&& f, Gradient&& g, HessianOp&& H, PrecondOp&& invB, const ASize size, Scalar* x) requires (isSize<ASize>) { return solveInPlace(makeOracle<Scalar>(size, std::forward<Function>(f), std::forward<Gradient>(g), std::forward<HessianOp>(H), std::forward<PrecondOp>(invB)), x); }

	/**
	 * @brief Solve in place using raw function, gradient, Hessian product and preconditioner operator functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param H Hessian-product functor.
	 * @param invB Preconditioner functor.
	 * @param x Initial guess and output vector for the solution.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient, CHessianOp<Scalar> HessianOp, CHessianOp<Scalar> PrecondOp>
	void solveInPlace(Function&& f, Gradient&& g, HessianOp&& H, PrecondOp&& invB, const std::span<Scalar> x) { return solveInPlace(makeOracle<Scalar>(Size(x.size()), std::forward<Function>(f), std::forward<Gradient>(g), std::forward<HessianOp>(H), std::forward<PrecondOp>(invB)), x); }
	
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
	void solveWithGuess(Function&& f, Gradient&& g, HessianOp&& H, PrecondOp&& invB, const Scalar* x0, const ASize size, Scalar* x) requires (isSize<ASize>) { return solveWithGuess(makeOracle<Scalar>(size, std::forward<Function>(f), std::forward<Gradient>(g), std::forward<HessianOp>(H), std::forward<PrecondOp>(invB)), x0, x); }

	/**
	 * @brief Solve with an initial guess using raw function, gradient, Hessian product and preconditioner operator functors.
	 * @param f Function functor.
	 * @param g Gradient functor.
	 * @param H Hessian-product functor.
	 * @param invB Preconditioner functor.
	 * @param x0 Initial guess.
	 * @param x Output solution vector.
	 */
	template<CFunction<Scalar> Function, CGradient<Scalar> Gradient, CHessianOp<Scalar> HessianOp, CHessianOp<Scalar> PrecondOp, std::ranges::input_range InitialGuess>
	void solveWithGuess(Function&& f, Gradient&& g, HessianOp&& H, PrecondOp&& invB, InitialGuess&& x0, const std::span<Scalar> x) { return solveWithGuess(makeOracle<Scalar>(Size(x.size()), std::forward<Function>(f), std::forward<Gradient>(g), std::forward<HessianOp>(H), std::forward<PrecondOp>(invB)), std::forward<InitialGuess>(x0), x); }

	// ===================================================================
	// MONITORING METHODS
	// ===================================================================
	
	Scalar getValue      () const { return m_fx;       } ///<  @brief Get last computed function value
	Scalar getError      () const { return m_residual; } ///<  @brief Get L2 norm of the last computed gradient
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

#define LNOT_DEFINE_SECOND_ORDER_SOLVER \
	using Base   = SecondOrderSolverBase<Self>; \
	using Size   = typename Base::Size; \
	using Scalar = typename Base::Scalar; \
	using Info   = typename Base::Info; \
	\
	template<typename ASize> static constexpr bool isSize = Base::template isSize<ASize>; \
	template<typename ABool> static constexpr bool isBool = Base::template isBool<ABool>; \

#define LNOT_SECOND_ORDER_SOLVER_ATTRIBUTE \
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

template<class Solver> concept CSecondOrderSolver = std::derived_from<Solver, SecondOrderSolverBase<Solver>>;

} // namespace LightOptim

#endif // LNOT_SECOND_ORDER_SOLVER_BASE_HPP
