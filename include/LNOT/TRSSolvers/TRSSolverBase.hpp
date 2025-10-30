#ifndef LNOT_TRS_SOLVER_BASE_HPP
#define LNOT_TRS_SOLVER_BASE_HPP

#include <cstdio>
#include <cmath>

#include <type_traits>
#include <limits>
#include <algorithm>

#include <LNOT/BasicLinalg/IdentityPreconditionerOp.hpp>

namespace LNOT
{
	
template<class TRSSolver> struct TRSSolverTraits;

template<class Derived>
class TRSSolverBase
{
	using DerivedTraits = TRSSolverTraits<Derived>;
public:
	using Scalar = typename DerivedTraits::Scalar;
	using Size   = typename DerivedTraits::Size;
	
	enum class Info {SUCCESS, FAILURE, BREAKDOWN};
	
	using IdOp = IdentityPreconditionerOp<Scalar>;
	
	template<typename HesOp>                  struct IsHessianOp   : std::bool_constant< std::is_invocable<HesOp, const Scalar*, Scalar*>::value > {}; ///<  @brief Trait to check if a type is a valid Hessian operator.
	template<typename HesOp, typename PrecOp> struct AreHessianOps : std::bool_constant< IsHessianOp<HesOp>::value and IsHessianOp<PrecOp>::value> {}; ///<  @brief Trait to check if two types are both valid Hessian operators.

	const Derived& derived() const { return static_cast<const Derived&>(*this); }
	      Derived& derived()       { return static_cast<      Derived&>(*this); }
	
	TRSSolverBase(const Size maxIt = 200000, const Scalar tol = std::numeric_limits<Scalar>::epsilon(), const Scalar tolTr = std::sqrt(std::numeric_limits<Scalar>::epsilon())) : m_maxIt(maxIt), m_tol(tol), m_tolTr(tolTr) {}
	~TRSSolverBase() { clearWorkSpace(); }
	      
	void clearWorkSpace() { derived().clearWorkSpace(); }
	
	template<typename Op> 
	void solve(const Op& H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) requires (IsHessianOp<Op>::value) { IdOp I(size); derived().solveImpl(H, I, g, size, delta, x); }
	
	template<typename HesOp, typename PrecOp> 
	void solve(const HesOp& H, const PrecOp& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) requires (AreHessianOps<HesOp, PrecOp>::value) { derived().solveImpl(H, invB, g, size, delta, x); }

	Scalar getError        () const { return derived().getErrorImpl();        }
	Scalar getSquaredError () const { return derived().getSquaredErrorImpl(); }
	
	Size   getMaxIt      () const { return m_maxIt; }
	Scalar getTol        () const { return m_tol;   }
	Scalar getTolTR      () const { return m_tolTr; }
	Size   getIterations () const { return m_nIt;   }
	Info   getInfo       () const { return m_info;  }
	
	Scalar getModelReduction() const { return m_modelReduction; }
	
	void setMaxIt (const Size    maxIt) { m_maxIt = maxIt; }
	void setTol   (const Scalar& tol)   { m_tol   = tol;   }
	void setTolTR (const Scalar& tolTr) { m_tolTr = tolTr; }
	
	void setOutput(std::FILE* out) { m_out = out; }
protected:
	Scalar getResidualThreshold()        const { return m_tol*std::max(Scalar(1), getError()); }              ///<  @brief Get the residual threshold for which the problem is considered solved.
	Scalar getSquaredResidualThreshold() const { return m_tol*m_tol*std::max(Scalar(1), getSquaredError()); } ///<  @brief Get the squared residual threshold for which the problem is considered solved. 

	Size   m_maxIt;
	Scalar m_tol;
	Scalar m_tolTr;
	Size   m_nIt;
	Info   m_info;
	
	Scalar m_modelReduction = 0;
	
	Size m_workCapacity = 0;
	
	std::FILE* m_out = nullptr;
};

template<class T> struct IsTRSSolver : std::bool_constant< std::is_base_of<TRSSolverBase<T>, T>::value > {};

template<class T> concept TRSSolver_concept = IsTRSSolver<T>::value;
	
} // namespace LNOT

#endif // LNOT_TRS_SOLVER_BASE_HPP
