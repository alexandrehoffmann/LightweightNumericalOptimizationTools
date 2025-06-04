#ifndef LNOT_LINEAR_SOLVER_BASE_HPP
#define LNOT_LINEAR_SOLVER_BASE_HPP

#include <cstdio>

#include <type_traits>
#include <limits>
#include <algorithm>

namespace LNOT
{

template<class LinearSolver> struct LinearSolverTraits;

template<class Derived>
class LinearSolverBase
{
	using DerivedTraits = LinearSolverTraits<Derived>;
public:
	using Scalar = typename DerivedTraits::Scalar;
	using Size   = typename DerivedTraits::Size;
	
	enum class Info {SUCCESS, FAILURE, NEGATIVE_CURVATURE, BREAKDOWN};
	
	template<typename HessianOp> struct IsHessianOp : std::bool_constant< std::is_invocable<HessianOp, const Scalar*, Scalar*>::value > {};

	const Derived& derived_cast() const { return static_cast<const Derived&>(*this); }
	      Derived& derived_cast()       { return static_cast<      Derived&>(*this); }
	
	LinearSolverBase(const Size maxIt = 200000, const Scalar tol = std::numeric_limits<Scalar>::epsilon()) : m_maxIt(maxIt), m_tol(tol) {}
	~LinearSolverBase() { clearWorkSpace(); }
	
	void clearWorkSpace() { derived_cast().clearWorkSpace(); }
	
	template<typename Op> 
	void solve(const Op& H, const Scalar* g, const Size size, Scalar* x) requires (IsHessianOp<Op>::value) { derived_cast().solve_impl(H, g, size, std::false_type{}, x); }
	
	template<typename Op> 
	void solveWithGuess(const Op& H, const Scalar* g, const Scalar* x0, const Size size, Scalar* x) requires (IsHessianOp<Op>::value) { std::copy(x0, x0 + size, x); derived_cast().solve_impl(H, g, size, std::true_type{}, x);  }
	
	template<typename Op, bool solveInPlace> 
	void solve_impl(const Op& H, const Scalar* g, const Size size, std::bool_constant<solveInPlace> bc, Scalar* x) requires (IsHessianOp<Op>::value) { derived_cast().solve_impl(H, g, size, bc, x); }
	
	Scalar getError        () const { return derived_cast().getError();        }
	Scalar getSquaredError () const { return derived_cast().getSquaredError(); }
	
	Size   getMaxIt      () const { return m_maxIt; }
	Scalar getTol        () const { return m_tol;   }
	Size   getIterations () const { return m_nIt;   }
	Info   getInfo       () const { return m_info;  }
	
	void setMaxIt (const Size    maxIt) { m_maxIt = maxIt; }
	void setTol   (const Scalar& tol)   { m_tol   = tol;   }
	
	void setOutput(std::FILE* out) { m_out = out; }
protected:
	Size   m_maxIt;
	Scalar m_tol;
	Size   m_nIt;
	Info   m_info;
	
	Size m_workCapacity = 0;
	
	std::FILE* m_out = nullptr;
};

template<class T> struct IsLinearSolver : std::bool_constant< std::is_base_of<LinearSolverBase<T>, T>::value > {};

template<class T> concept LinearSolver_concept = IsLinearSolver<T>::value;

} // namespace LNOT

#endif // LNOT_LINEAR_SOLVER_BASE_HPP
