#ifndef LNOT_TRS_SOLVER_BASE_HPP
#define LNOT_TRS_SOLVER_BASE_HPP

#include <cstdio>
#include <cmath>

#include <type_traits>
#include <algorithm>

#include <LNOT/BasicLinalg/IdentityPreconditionerOp.hpp>
#include <LNOT/FloatingPoint/NumTraits.hpp>
#include <LNOT/CRTPBase.hpp>
#include <LNOT/misc/AdlMath.hpp>

#include <BIC/Core.hpp>

namespace LNOT
{
	
template<class TRSSolver> struct TRSSolverTraits;

template<class Derived>
class TRSSolverBase : public CRTPBase<Derived>
{
	using DerivedTraits = TRSSolverTraits<Derived>;
public:
	using Scalar = typename DerivedTraits::Scalar; ///<  @brief The scalar type used in computations (e.g., float, double)
	using Size   = typename DerivedTraits::Size;   ///<  @brief The size type used for indexing and loop counters
	
	enum class Info {SUCCESS, FAILURE, BREAKDOWN}; ///<  @brief Enumeration indicating solver termination status.
	
	using IdOp = IdentityPreconditionerOp<Scalar>;
	
	template<typename HesOp>                  struct IsHessianOp   : BIC::Fixed<bool, std::is_invocable<HesOp, const Scalar*, Scalar*>::value > {}; ///<  @brief Trait to check if a type is a valid Hessian operator.
	template<typename HesOp, typename PrecOp> struct AreHessianOps : BIC::Fixed<bool, IsHessianOp<HesOp>::value and IsHessianOp<PrecOp>::value> {}; ///<  @brief Trait to check if two types are both valid Hessian operators.
	template<typename ASize>                  struct IsSize        : BIC::Fixed<bool, std::is_same<Size, BIC::Mutable<ASize>>::value > {};          ///<  @brief Trait to check if a type is either a `Size` or a `BIC::Fixed<Size, VALUE>`.
	
	static constexpr Scalar defaultEps = NumTraits<Scalar>::epsilon; ///<  @brief Default value for relative and absolute tolerance of the solver.  
	
	TRSSolverBase(const Size maxIt = 200000, const Scalar relTol = defaultEps, const Scalar relTolTr = AdlMath::sqrt(defaultEps), const Scalar absTol = defaultEps, const Scalar absTolTr = AdlMath::sqrt(defaultEps)) : m_maxIt(maxIt), m_relTol(relTol), m_relTolTr(relTolTr), m_absTol(absTol), m_absTolTr(absTolTr) {}
	
	template<typename Op, typename ASize> 
	Scalar solve(const Op& H, const Scalar* g, const ASize size, const Scalar& delta, Scalar* x) requires (IsHessianOp<Op>::value and IsSize<ASize>::value) { IdOp I(size); return this->derived().solveImpl(H, I, g, size, delta, x); }
	
	template<typename HesOp, typename PrecOp, typename ASize> 
	Scalar solve(const HesOp& H, const PrecOp& invB, const Scalar* g, const ASize size, const Scalar& delta, Scalar* x) requires (AreHessianOps<HesOp, PrecOp>::value and IsSize<ASize>::value) { return this->derived().solveImpl(H, invB, g, size, delta, x); }

	Scalar getError        () const { return this->derived().getErrorImpl();        }
	Scalar getSquaredError () const { return this->derived().getSquaredErrorImpl(); }
	
	Size   getMaxIt      () const { return m_maxIt;    }
	Scalar getRelTol     () const { return m_relTol;   }
	Scalar getRelTolTR   () const { return m_relTolTr; }
	Scalar getAbsTol     () const { return m_absTol;   }
	Scalar getAbsTolTR   () const { return m_absTolTr; }
	Size   getIterations () const { return m_nIt;      }
	Info   getInfo       () const { return m_info;     }
	
	Scalar getModelReduction() const { return m_modelReduction; }
	
	void setMaxIt    (const Size    maxIt) { m_maxIt    = maxIt; }
	void setRelTol   (const Scalar& tol)   { m_relTol   = tol;   }
	void setRelTolTR (const Scalar& tolTr) { m_relTolTr = tolTr; }
	void setAbsTol   (const Scalar& tol)   { m_absTol   = tol;   }
	void setAbsTolTR (const Scalar& tolTr) { m_absTolTr = tolTr; }
	
	void setTol   (const Scalar& tol)   { setRelTol(tol);     setAbsTol(tol);     } ///<  @brief Set both relative and absolute tolerance.
	void setTolTR (const Scalar& tolTr) { setRelTolTR(tolTr); setAbsTolTR(tolTr); } ///<  @brief Set both relative and absolute tolerance.
	
	void setOutput(std::FILE* out) { m_out = out; }
protected:
	Size   m_maxIt;
	Scalar m_relTol;
	Scalar m_relTolTr;
	Scalar m_absTol;
	Scalar m_absTolTr;
	Size   m_nIt;
	Info   m_info;
	
	Scalar m_modelReduction = Scalar{};
	
	Size m_workCapacity = Size{};
	
	std::FILE* m_out = nullptr;
};

#define LNOT_DEFINE_TRS_SOLVER \
	using Base   = TRSSolverBase<Self>; \
	using Size   = typename Base::Size; \
	using Scalar = typename Base::Scalar; \
	using Info   = typename Base::Info; \
	\
	template<typename HesOp, typename PrecOp> using AreHessianOps = typename Base::template AreHessianOps<HesOp,PrecOp>; \
	template<typename ASize>                  using IsSize        = typename Base::template IsSize<ASize>; \

#define LNOT_TRS_SOLVER_ATTRIBUTE \
	using Base::m_maxIt; \
	using Base::m_relTol; \
	using Base::m_relTolTr; \
	using Base::m_absTol; \
	using Base::m_absTolTr; \
	using Base::m_nIt; \
	using Base::m_info; \
	using Base::m_modelReduction; \
	using Base::m_workCapacity; \
	using Base::m_out; \

template<class T> struct IsTRSSolver : BIC::Fixed<bool, std::is_base_of<TRSSolverBase<T>, T>::value > {};

template<class T> concept CTRSSolver = IsTRSSolver<T>::value;
	
} // namespace LNOT

#endif // LNOT_TRS_SOLVER_BASE_HPP
