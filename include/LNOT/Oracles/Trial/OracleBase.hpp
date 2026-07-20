#ifndef LNOT_ORACLE_BASE_HPP
#define LNOT_ORACLE_BASE_HPP

#include <LNOT/CRTPBase.hpp>
#include <LNOT/Oracles/concepts.hpp>
#include <LNOT/Oracles/helpers.hpp>

#include <BIC/Core.hpp>

namespace LNOT
{

template<class Derived> struct OracleTraits;

template<class Derived>
class OracleBase : public CRTPBase<Derived>
{
	using DerivedTraits = OracleTraits<Derived>;
	using CRTP		  = CRTPBase<Derived>;
public:
	using Size	  = typename DerivedTraits::Size;   ///<  @brief Type representing the number of dimensions (variables) of the function.
	using Scalar	= typename DerivedTraits::Scalar; ///<  @brief Scalar type used in evaluations.
	using HessianOp = typename detail::HessianOpForTrait<DerivedTraits>::Type;
	using PrecondOp = typename detail::PrecondOpForTrait<DerivedTraits>::Type;
	
	static constexpr bool hasGradient  = CFirstOrderInterface<Derived, Scalar>;
	static constexpr bool hasHessianOp = CSecondOrderInterface<Derived, Scalar, HessianOp>;
	static constexpr bool hasPrecondOp = CPrecondFirstOrderInterface<Derived, Scalar, PrecondOp> or CPrecondhSecondOrderInterface<Derived, Scalar, HessianOp, PrecondOp>;
	
	/**
	 * @brief Get the number of dimensions (variables) of the function.
	 * 
	 * Delegates to `Derived::getNDimsImpl()`.
	 * 
	 * @return Number of variables (dimension of the domain).
	 */
	constexpr Size getNDims() const { return CRTP::derived().getNDimsImpl(); }
	
	constexpr HessianOp makeHessianOp() requires(hasHessianOp) { return CRTP::derived().makeHessianOpImpl(); }
	
	constexpr PrecondOp makePrecondOp() requires(hasGradient);
	
	constexpr Scalar eval(const Scalar* x) { ++m_nFuncEval; return CRTP::derived().getValue(x); }
	
	constexpr Scalar eval(const Scalar* x, Scalar* g) requires(hasGradient);
	
	constexpr Scalar eval(const Scalar* x, Scalar* g, PrecondOp& invBop) requires(hasGradient);
	
	constexpr Scalar eval(const Scalar* x, Scalar* g, HessianOp& Hop) requires(hasHessianOp);
	
	constexpr Scalar eval(const Scalar* x, Scalar* g, HessianOp& Hop, PrecondOp& invBop) requires(hasHessianOp);
	
	constexpr Size getNFuncEval() const { return m_nFuncEval; }
	
	constexpr Size getNGradEval() const requires(hasGradient) { return m_nGradEval; }
	
	constexpr Size getNHessianOpEvalEval() const requires(hasHessianOp) { return m_nHessianOpEval; }
	
	constexpr Size getNPrecondOpEval() const requires(hasPrecondOp) { return m_nPrecondOpEval; }
private:
	Size m_nFuncEval      = {};
	Size m_nGradEval      = {};
	Size m_nHessianOpEval = {};
	Size m_nPrecondOpEval = {};
};

#define LNOT_DEFINE_ORACLE \
	using Base      = LNOT::OracleBase<Self>; \
	using Size      = typename Base::Size; \
	using Scalar    = typename Base::Scalar; \
	using HessianOp = typename Base::HessianOp; \
	using PrecondOp = typename Base::PrecondOp; \
	\
	using Base::hasGradient; \
	using Base::hasHessianOp; \
	using Base::hasPrecondOp; \

template<class Oracle> concept COracle            = std::derived_from<std::remove_cvref_t<Oracle>, OracleBase<std::remove_cvref_t<Oracle>>>;
template<class Oracle> concept CFirstOrderOracle  = COracle<Oracle> and std::remove_cvref_t<Oracle>::hasGradient;
template<class Oracle> concept CSecondOrderOracle = COracle<Oracle> and std::remove_cvref_t<Oracle>::hasGradient and std::remove_cvref_t<Oracle>::hasHessianOp;

} // namespace LNOT

#include <LNOT/Oracles/OracleBase_impl.hpp>

#endif // LNOT_ORACLE_BASE_HPP

