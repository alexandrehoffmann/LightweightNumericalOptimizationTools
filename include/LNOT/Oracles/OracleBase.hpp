#ifndef LNOT_ORACLE_BASE_HPP
#define LNOT_ORACLE_BASE_HPP

#include <LNOT/Traits.hpp>
#include <LNOT/CRTPBase.hpp>
#include <LNOT/Oracles/concepts.hpp>
#include <LNOT/Oracles/HessianOpAdaptor.hpp>
#include <LNOT/Oracles/PrecondOpAdaptor.hpp>
#include <LNOT/BasicLinalg/IdentityOp.hpp>

#include <BIC/Core.hpp>

#include <type_traits> // std::bool_constant and std::is_base_of
#include <algorithm>   // std::copy

namespace LNOT
{

template<class Derived> struct OracleTraits;

/**
 * @brief Base class for an Oracle using the CRTP (Curiously Recurring Template Pattern).
 * 
 * An Oracle provides basic information about a real-valued function \f$f : \mathbb{R}^d \to \mathbb{R}\f$ at a given point \f$x \in \mathbb{R}^d\f$:
 * - getValue() returns \f$f(x)\f$
 * - getGradient() (optional) fills g with the values of \f$\nabla f(x)\f$
 * - getHessianProd() (optional) fills Hd with the values of \f$\nabla^2 f(x) d\f$
 * - applyPrecond() (optional) apply a preconditioner. If the derived class doesn't implement applyPrecond(), the Identity preconditioner will be used
 * 
 * LNOT's Oracle work in the following way. 
 * First, a \f$x\in\mathbb{R}^d\f$ is provided via the method setCurrentPoint(). 
 * The Oracle <em>does not</em> take ownership of the array x. However, it may use it to do some pre-computations needed to evaluation 
 * function's value, gradient and/or Hessian product. 
 * \em e.g.  for a frequency domaine Full Waveform Inversion problem, setCurrentPoint() would basically build and factorize the problem mass and stiffness matrices.
 * 
 * @tparam Derived The derived class that inherits from OracleBase. 
 * Derived must implement:
 * - getNDimsImpl()
 * - setCurrentPointImpl()
 * - getValueImpl()
 * It may also implement:
 * - getGradientImpl()
 * - getHessianProdImpl() and initHessianProd()
 * - applyPrecondImpl() and initApplyPrecond()
 * An OracleTraits must be specialized for Derived that defines:
 * - Scalar, the Scalar type used by the Oracle, typically float or double
 * - Size, the Size type, typically int, unsigned int, size_t
 */
template<class Derived>
class OracleBase : public CRTPBase<Derived>
{
	using DerivedTraits = OracleTraits<Derived>;
	using CRTP          = CRTPBase<Derived>;
public:
	using Size   = typename DerivedTraits::Size;   ///<  @brief Type representing the number of dimensions (variables) of the function.
	using Scalar = typename DerivedTraits::Scalar; ///<  @brief Scalar type used in evaluations.
	/**
	 * @brief Get the number of dimensions (variables) of the function.
	 * 
	 * Delegates to `Derived::getNDimsImpl()`.
	 * 
	 * @return Number of variables (dimension of the domain).
	 */
	constexpr std::convertible_to<Size> auto getNDims() const requires(COracleInterface<Derived, Scalar, Size>) { return CRTP::derived().getNDimsImpl(); }
	
	/**
	 * @brief Set the point at which the function will be evaluated.
	 * 
	 * This method does not take ownership of the pointer.
	 * Delegates to `Derived::setCurrentPointImpl(x)`.
	 * 
	 * @param x Pointer to an array representing the current point \f$x \in \mathbb{R}^d\f$.
	 */
	constexpr void setCurrentPoint(const Scalar* x) requires(COracleInterface<Derived, Scalar, Size>) { CRTP::derived().setCurrentPointImpl(x); }
	
	/**
	 * @brief True if the current point is feasible
	 * 
	 * Delegates to `Derived::isFeasibleImpl(x)`.
	 * 
	 * @return true if the current point is feasible.
	 */
	 constexpr bool isFeasible() const requires(COracleInterface<Derived, Scalar, Size>) { return CRTP::derived().isFeasibleImpl(); }
	
	/**
	 * @brief Get the function value \f$f(x)\f$ at the current point.
	 * 
	 * Delegates to `Derived::getValueImpl()`.
	 * 
	 * @return The function value at the current point.
	 */
	constexpr std::convertible_to<Scalar> auto getValue() const requires(COracleInterface<Derived, Scalar, Size>) { return CRTP::derived().getValueImpl(); }
	
	/**
	 * @brief Get the gradient \f$\nabla f(x)\f$ at the current point.
	 * 
	 * Enabled only if `hasGradient == true`.
	 * 
	 * Delegates to `Derived::getGradientImpl(g)`.
	 * 
	 * @param g Pointer to an array to store the gradient values.
	 */
	constexpr void getGradient(Scalar* g) const requires(CFirstOrderOracleInterface<Derived, Scalar, Size>) { CRTP::derived().getGradientImpl(g); }
	
	constexpr HessianOpAdaptor<Derived> getHessianOp() requires(CSecondOrderOracleInterface<Derived, Scalar, Size>) { return HessianOpAdaptor(CRTP::derived()); }
	
	constexpr PrecondOpAdaptor<Derived> getPrecondOp() requires(CPrecondOracleInterface<Derived, Scalar>) { return PrecondOpAdaptor(CRTP::derived()); }
	
	constexpr IdentityOp getPrecondOp() const requires(not CPrecondOracleInterface<Derived, Scalar>) { return IdentityOp(getNDims()); }
};

#define LNOT_DEFINE_ORACLE \
	using Base   = LNOT::OracleBase<Self>; \
	using Size   = typename Base::Size; \
	using Scalar = typename Base::Scalar; \

template<class Oracle> concept COracle            = std::derived_from<std::remove_cvref_t<Oracle>, OracleBase<std::remove_cvref_t<Oracle>>>;
template<class Oracle> concept CFirstOrderOracle  = COracle<Oracle> and CFirstOrderOracleInterface<Oracle, ScalarFor<Oracle>, SizeFor<Oracle>>;
template<class Oracle> concept CSecondOrderOracle = COracle<Oracle> and CSecondOrderOracleInterface<Oracle, ScalarFor<Oracle>, SizeFor<Oracle>>;

} // namespace LNOT

#endif // LNOT_ORACLE_BASE_HPP

