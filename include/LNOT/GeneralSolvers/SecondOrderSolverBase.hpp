#ifndef LNOT_SECOND_ORDER_MINIMIZER_BASE_HPP
#define LNOT_SECOND_ORDER_MINIMIZER_BASE_HPP

#include <LNOT/Oracle/OracleWrapper.hpp>

#include <cstdio>
#include <cmath>

#include <vector>
#include <algorithm>

namespace LNOT
{

template<class Derived> struct SecondOrderSolverTraits;

template<class Derived>
class SecondOrderSolverBase
{
	using DerivedTraits = SecondOrderSolverTraits<Derived>;
public:
	using Scalar = typename DerivedTraits::Scalar;
	using Size   = typename DerivedTraits::Size;
	
	enum class Info {SUCCESS, FAILURE, BREAKDOWN};
	
	const Derived& derived_cast() const { return static_cast<const Derived&>(*this); }
	      Derived& derived_cast()       { return static_cast<      Derived&>(*this); }
	
	template<typename Function>  struct IsFunction  : std::bool_constant< std::is_invocable<Function, const Scalar*>::value > {};
	template<typename Gradient>  struct IsGradient  : std::bool_constant< std::is_invocable<Gradient, const Scalar*, Scalar*>::value > {};
	template<typename HessianOp> struct IsHessianOp : std::bool_constant< std::is_invocable<HessianOp, const Scalar*, const Scalar*, Scalar*>::value > {};
	template<typename Function, typename Gradient, typename HessianOp> struct IsProgram : std::bool_constant< IsFunction<Function>::value and IsGradient<Gradient>::value and IsHessianOp<HessianOp>::value > {};
	//////////////////////////////////////////////////////////////////////
	SecondOrderSolverBase(const Size maxIt = 200000, const Scalar tol = std::numeric_limits<Scalar>::epsilon()) : m_maxIt(maxIt), m_tol(tol) {}
	~SecondOrderSolverBase() { clearWorkSpace(); }
	
	void clearWorkSpace() { derived_cast().clearWorkSpace(); }
	//////////////////////////////////////////////////////////////////////
	template<SecondOrderOracle_concept Oracle>
	void solve(Oracle& oracle, Scalar* x) { solve_impl(oracle, std::false_type{}, x); }
	
	template<SecondOrderOracle_concept Oracle>
	void solveWithGuess(Oracle& oracle, const Scalar* x0, Scalar* x) { std::copy(x0, x0 + oracle.getNDims(), x); solve_impl(oracle, std::true_type{}, x); }
	//////////////////////////////////////////////////////////////////////
	template<typename Function, typename Gradient, typename HessianOp>
	void solve(Function f, Gradient g, HessianOp H, const Size size, Scalar* x) requires (IsProgram<Function,Gradient,HessianOp>::value) { solve_impl(f, g, H, size, std::false_type{}, x); }
	
	template<typename Function, typename Gradient, typename HessianOp>
	void solveWithGuess(Function f, Gradient g, HessianOp H, const Scalar* x0, const Size size, Scalar* x) requires (IsProgram<Function,Gradient,HessianOp>::value) { std::copy(x0, x0 + size, x); solve_impl(f, g, H, size, std::true_type{}, x);  }
	//////////////////////////////////////////////////////////////////////
	template<typename Function, typename Gradient, typename HessianOp, bool solveInPlace> 
	void solve_impl(Function f, Gradient g, HessianOp H, const Size size, std::bool_constant<solveInPlace> bc, Scalar* x) requires (IsProgram<Function,Gradient,HessianOp>::value) { OracleWrapper<Scalar,Function,Gradient,HessianOp> oracle(size, f, g, H); solve_impl(oracle, bc, x); }
	
	template<SecondOrderOracle_concept Oracle, bool solveInPlace> 
	void solve_impl(Oracle& oracle, std::bool_constant<solveInPlace> bc, Scalar* x) { derived_cast().solve_impl(oracle, bc, x); }
	//////////////////////////////////////////////////////////////////////
	Scalar getValue        () const { return m_fx; }
	Scalar getError        () const { return std::sqrt(m_squaredNormGrad); }
	Scalar getSquaredError () const { return m_squaredNormGrad; }
	
	Size   getMaxIt      () const { return m_maxIt; }
	Scalar getTol        () const { return m_tol;   }
	Size   getIterations () const { return m_nIt;   }
	Info   getInfo       () const { return m_info;  }
	
	Size getInnerIterations(const Size it) const { return m_innerIts[it]; }
	
	void setMaxIt (const Size    maxIt) { m_maxIt = maxIt; }
	void setTol   (const Scalar& tol)   { m_tol   = tol;   }
	
	void setOutput(std::FILE* out) { m_out = out; }
protected:
	Size   m_maxIt;
	Scalar m_tol;
	Size   m_nIt;
	Info   m_info;
	Scalar m_squaredNormGrad;
	Scalar m_fx;
	
	Size m_workCapacity = 0;
	
	std::vector< Size > m_innerIts;
	
	std::FILE* m_out = nullptr;
};

template<class T> struct IsSecondOrderSolver : std::bool_constant< std::is_base_of<SecondOrderSolverBase<T>, T>::value > {};

template<class T> concept IsSecondOrderSolver_concept = IsSecondOrderSolver<T>::value;

} // namespace LightOptim

#endif // LNOT_SECOND_ORDER_MINIMIZER_BASE_HPP
