#ifndef LNOT_ORACLE_BASE_HPP
#define LNOT_ORACLE_BASE_HPP

#include <LNOT/CRTPBase.hpp>

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
 * - getHessianProdImpl()
 * - applyPrecondImpl()
 * An OracleTraits must be specialized for Derived that defines:
 * - Scalar, the Scalar type used by the Oracle, typically float or double
 * - Size, the Size type, typically int, unsigned int, size_t
 * - hasGradient, a boolean constexpr true, if Derived implements getGradient()
 * - hasHessianProd, a boolean constexpr true, if Derived implements getHessianProd()
 * - hasApplyPrecond, a boolean constexpr true, if Derived implements applyPrecond()
 */
template<class Derived>
class OracleBase : public CRTPBase<Derived>
{
	using DerivedTraits = OracleTraits<Derived>;
	using CRTP          = CRTPBase<Derived>;
public:
	using Size   = typename DerivedTraits::Size;   ///<  @brief Type representing the number of dimensions (variables) of the function.
	using Scalar = typename DerivedTraits::Scalar; ///<  @brief Scalar type used in evaluations.

	static constexpr bool hasGradient     = DerivedTraits::hasGradient;
	static constexpr bool hasHessianProd  = DerivedTraits::hasHessianProd;
	static constexpr bool hasApplyPrecond = DerivedTraits::hasApplyPrecond;
	
	/**
	 * @brief Get the number of dimensions (variables) of the function.
	 * 
	 * Delegates to `Derived::getNDimsImpl()`.
	 * 
	 * @return Number of variables (dimension of the domain).
	 */
	Size getNDims() const { return CRTP::derived().getNDimsImpl(); }
	
	/**
	 * @brief Set the point at which the function will be evaluated.
	 * 
	 * This method does not take ownership of the pointer.
	 * Delegates to `Derived::setCurrentPointImpl(x)`.
	 * 
	 * @param x Pointer to an array representing the current point \f$x \in \mathbb{R}^d\f$.
	 */
	void setCurrentPoint(const Scalar* x) { CRTP::derived().setCurrentPointImpl(x); }
	
	/**
	 * @brief True if the current point is feasible
	 * 
	 * Delegates to `Derived::isFeasibleImpl(x)`.
	 * 
	 * @return true if the current point is feasible.
	 */
	 bool isFeasible() const { return CRTP::derived().isFeasibleImpl(); }
	
	/**
	 * @brief Get the function value \f$f(x)\f$ at the current point.
	 * 
	 * Delegates to `Derived::getValueImpl()`.
	 * 
	 * @return The function value at the current point.
	 */
	Scalar getValue() const { return CRTP::derived().getValueImpl(); }
	
	/**
	 * @brief Get the gradient \f$\nabla f(x)\f$ at the current point.
	 * 
	 * Enabled only if `hasGradient == true`.
	 * 
	 * Delegates to `Derived::getGradientImpl(g)`.
	 * 
	 * @param g Pointer to an array to store the gradient values.
	 */
	void getGradient(Scalar* g) const requires (hasGradient) { CRTP::derived().getGradientImpl(g); }
	
	/**
	 * @brief Compute the Hessian-vector product \f$\nabla^2 f(x)d\f$.
	 * 
	 * Enabled only if `hasHessianProd == true`.
	 * 
	 * Delegates to `Derived::getHessianProdImpl(d, Hd)`.
	 * 
	 * @param d  Input direction vector.
	 * @param Hd Output array for the result of the Hessian-vector product.
	 */
	void getHessianProd(const Scalar* d, Scalar* Hd) const requires (hasHessianProd) { CRTP::derived().getHessianProdImpl(d, Hd); }
	
	/**
	 * @brief Approximatively solves \f$B x = y\f$ where \f$B\approx\nabla^2 f(x)\f$.
	 * 
	 * Enabled only if `hasHessianProd == true`.
	 * 
	 * Delegates to `Derived::hasApplyPrecondImpl(d, invBd)`.
	 * 
	 * @param d  Input direction vector.
	 * @param invBd Output array for the result.
	 */
	void applyPrecond(const Scalar* d, Scalar* invBd) const requires (hasApplyPrecond) { CRTP::derived().applyPrecondImpl(d, invBd); }
	
	/**
	 * @brief dafault implementation for applyPrecondImpl. Copy d into invBd.
	 */
	void applyPrecond(const Scalar* d, Scalar* invBd) const requires (not hasApplyPrecond) { std::copy(d, d + getNDims(), invBd); }
};

#define LNOT_DEFINE_ORACLE \
	using Base   = LNOT::OracleBase<Self>; \
	using Size   = typename Base::Size; \
	using Scalar = typename Base::Scalar; \
	\
	using Base::hasGradient; \
	using Base::hasHessianProd; \
	using Base::hasApplyPrecond; \

template<class Oracle> concept COracle            = std::derived_from<Oracle, OracleBase<Oracle>>;
template<class Oracle> concept CFirstOrderOracle  = COracle<Oracle> and Oracle::hasGradient;
template<class Oracle> concept CSecondOrderOracle = COracle<Oracle> and Oracle::hasGradient and Oracle::hasHessianProd;

} // namespace LNOT

#endif // LNOT_ORACLE_BASE_HPP

