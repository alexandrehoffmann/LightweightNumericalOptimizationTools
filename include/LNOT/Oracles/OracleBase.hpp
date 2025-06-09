#ifndef LNOT_ORACLE_BASE_HPP
#define LNOT_ORACLE_BASE_HPP

#include <type_traits> // std::bool_constant and std::is_base_of

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
 * 
 * LNOT's Oracle work in the following way. 
 * First, a \f$x\in\mathbb{R}^d\f$ is provided via the method setCurrentPoint(). 
 * The Oracle <em>does not</em> take ownership of the array x. However, it may use it to do some pre-computations needed to evaluation 
 * function's value, gradient and/or Hessian product. 
 * \em e.g.  for a frequency domaine Full Waveform Inversion problem, setCurrentPoint() would basically build and factorize the problem mass and stiffness matrices.
 * 
 * @tparam Derived The derived class that inherits from OracleBase. 
 * Derived must implement:
 * - getNDims()
 * - setCurrentPoint()
 * - getValue()
 * It may also implement:
 * - getGradient()
 * - getHessianProd()
 * An OracleTraits must be specialized for Derived that defines:
 * - Scalar, the Scalar type used by the Oracle, typically float or double
 * - Size, the Size type, typically int, unsigned int, size_t
 * - hasGradient, a boolean constexpr true, if Derived implements getGradient()
 * - hasHessianProd, a boolean constexpr true, if Derived implements getHessianProd()
 */
template<class Derived>
class OracleBase
{
	using DerivedTraits = OracleTraits<Derived>;
public:
	using Size   = typename DerivedTraits::Size;   ///<  @brief Type representing the number of dimensions (variables) of the function.
	using Scalar = typename DerivedTraits::Scalar; ///<  @brief Scalar type used in evaluations.

	static constexpr bool hasGradient    = DerivedTraits::hasGradient;
	static constexpr bool hasHessianProd = DerivedTraits::hasHessianProd;

	const Derived& derived_cast() const { return static_cast<const Derived&>(*this); }
	      Derived& derived_cast()       { return static_cast<      Derived&>(*this); }
	/**
	 * @brief Get the number of dimensions (variables) of the function.
	 * 
	 * Delegates to `Derived::getNDims()`.
	 * 
	 * @return Number of variables (dimension of the domain).
	 */
	Size getNDims() const { return derived_cast().getNDims(); }
	/**
	 * @brief Set the point at which the function will be evaluated.
	 * 
	 * This method does not take ownership of the pointer.
	 * Delegates to `Derived::setCurrentPoint(x)`.
	 * 
	 * @param x Pointer to an array representing the current point \f$x \in \mathbb{R}^d\f$.
	 */
	void setCurrentPoint(const Scalar* x) { derived_cast().setCurrentPoint(x); }
	/**
	 * @brief Get the function value \f$f(x)\f$ at the current point.
	 * 
	 * Delegates to `Derived::getValue()`.
	 * 
	 * @return The function value at the current point.
	 */
	Scalar getValue() const { return derived_cast().getValue(); }
	/**
	 * @brief Get the gradient \f$\nabla f(x)\f$ at the current point.
	 * 
	 * Enabled only if `hasGradient == true`.
	 * 
	 * Delegates to `Derived::getGradient(g)`.
	 * 
	 * @param g Pointer to an array to store the gradient values.
	 */
	void getGradient(Scalar* g) const requires (hasGradient) { derived_cast().getGradient(g); }
	/**
	 * @brief Compute the Hessian-vector product \f$\nabla^2 f(x)d\f$.
	 * 
	 * Enabled only if `hasHessianProd == true`.
	 * 
	 * Delegates to `Derived::getHessianProd(d, Hd)`.
	 * 
	 * @param d  Input direction vector.
	 * @param Hd Output array for the result of the Hessian-vector product.
	 */
	void getHessianProd(const Scalar* d, Scalar* Hd) const requires (hasHessianProd) { derived_cast().getHessianProd(d, Hd); }
};

template<class T> struct IsOracle : std::bool_constant< std::is_base_of<OracleBase<T>, T>::value > {};

template<class T> concept Oracle_concept            = IsOracle<T>::value;
template<class T> concept FirstOrderOracle_concept  = IsOracle<T>::value and T::hasGradient;
template<class T> concept SecondOrderOracle_concept = IsOracle<T>::value and T::hasGradient and T::hasHessianProd;

} // namespace LNOT

#endif // LNOT_ORACLE_BASE_HPP

