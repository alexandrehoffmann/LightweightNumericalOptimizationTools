#ifndef LNOT_FIRST_ORDER_MINIMIZER_BASE_HPP
#define LNOT_FIRST_ORDER_MINIMIZER_BASE_HPP

#include <LNOT/Oracle/OracleWrapper.hpp>

namespace LNOT
{

template<class Derived> struct FirstOrderSolverTraits;

template<class Derived>
class FirstOrderSolverBase
{
	using DerivedTraits = FirstOrderSolverTraits<Derived>;
public:
	using Scalar = typename DerivedTraits::Scalar;
	using Size   = typename DerivedTraits::Size;
	
	enum class Info {SUCCESS, FAILURE, BREAKDOWN};
	
	const Derived& derived_cast() const { return static_cast<const Derived&>(*this); }
	      Derived& derived_cast()       { return static_cast<      Derived&>(*this); }
	
	template<typename Function>  struct IsFunction  : std::bool_constant< std::is_invocable<Function, const Scalar*>::value > {};
	template<typename Gradient>  struct IsGradient  : std::bool_constant< std::is_invocable<Gradient, const Scalar*, Scalar*>::value > {};
	template<typename Function, typename Gradient> struct IsProgram : std::bool_constant< IsFunction<Function>::value and IsGradient<Gradient>::value > {};
	//////////////////////////////////////////////////////////////////////
	FirstOrderSolverBase(const Size maxIt = 200000, const Scalar tol = std::numeric_limits<Scalar>::epsilon()) : m_maxIt(maxIt), m_tol(tol) {}
	~FirstOrderSolverBase() { clearWorkSpace(); }
	
	void clearWorkSpace() { derived_cast().clearWorkSpace(); }
	//////////////////////////////////////////////////////////////////////
	template<FirstOrderOracle_concept Oracle>
	void solve(Oracle& oracle, Scalar* x) { solve_impl(oracle, std::false_type{}, x); }
	
	template<FirstOrderOracle_concept Oracle>
	void solveWithGuess(Oracle& oracle, const Scalar* x0, Scalar* x) { std::copy(x0, x0 + oracle.getNDims(), x); solve_impl(oracle, std::true_type{}, x); }
	//////////////////////////////////////////////////////////////////////
	template<typename Function, typename Gradient>
	void solve(Function f, Gradient g, const Size size, Scalar* x) requires (IsProgram<Function,Gradient>::value) { solve_impl(f, g, size, std::false_type{}, x); }
	
	template<typename Function, typename Gradient>
	void solveWithGuess(Function f, Gradient g, const Scalar* x0, const Size size, Scalar* x) requires (IsProgram<Function,Gradient>::value) { std::copy(x0, x0 + size, x); solve_impl(f, g, size, std::true_type{}, x);  }
	//////////////////////////////////////////////////////////////////////
	template<typename Function, typename Gradient, bool solveInPlace> 
	void solve_impl(Function f, Gradient g, const Size size, std::bool_constant<solveInPlace> bc, Scalar* x) requires (IsProgram<Function,Gradient>::value) { OracleWrapper<Scalar,Function,Gradient> oracle(size, f, g); solve_impl(oracle, bc, x); }
	
	template<FirstOrderOracle_concept Oracle, bool solveInPlace> 
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

template<class T> struct IsFirstOrderSolver : std::bool_constant< std::is_base_of<FirstOrderSolverBase<T>, T>::value > {};

template<class T> concept IsFirstOrderSolver_concept = IsFirstOrderSolver<T>::value;

} // namespace LightOptim

#endif // LNOT_FIRST_ORDER_MINIMIZER_BASE_HPP
