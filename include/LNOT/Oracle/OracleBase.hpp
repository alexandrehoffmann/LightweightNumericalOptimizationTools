#ifndef LNOT_ORACLE_BASE_HPP
#define LNOT_ORACLE_BASE_HPP

#include <type_traits> // std::bool_constant and std::is_base_of

namespace LNOT
{

template<class Derived> struct OracleTraits;

/**
 * An Oracle provides basic informations about a real-valued function \f$f : \mathbb{R}^d\to\mathbb{R}$\f at a given point \f$x\in\mathbb{R}^d\f$: 
 * - getValue() returns \f$f(x)\f$
 * - getGradient() (optional) fills g with the values of \f$\nabla f(x)\f$
 * - getHessianProd() (optional) fills Hd with the values of \f$\nabla^2 f(x) d\f$
 * 
 * LNOT's Oracle work in the following way. 
 * First, a \f$x\in\mathbb{R}^d\f$ is provided via the method setCurrentPoint(x). 
 * The Oracle <em>does not</em> take ownership of the array x. However, it may use it to do some pre-computations needed to evaluation 
 * function's value, gradient and/or Hessian product. 
 * \em e.g.  for a frequency domaine Full Waveform Inversion problem, setCurrentPoint(x) would basically build and factorize the problem mass and stiffness matrices.
 */
template<class Derived>
class OracleBase
{
	using DerivedTraits = OracleTraits<Derived>;
public:
	using Size   = typename DerivedTraits::Size;
	using Scalar = typename DerivedTraits::Scalar;

	static constexpr bool hasGradient    = DerivedTraits::hasGradient;
	static constexpr bool hasHessianProd = DerivedTraits::hasHessianProd;

	const Derived& derived_cast() const { return static_cast<const Derived&>(*this); }
	      Derived& derived_cast()       { return static_cast<      Derived&>(*this); }
	      
	Size getNDims() const { return derived_cast().getNDims(); }

	void setCurrentPoint(const Scalar* x) { derived_cast().setCurrentPoint(x); }

	Scalar getValue() const { return derived_cast().getValue(); }

	void getGradient(Scalar* g) const requires (hasGradient) { derived_cast().getGradient(g); }

	void getHessianProd(const Scalar* d, Scalar* Hd) const requires (hasHessianProd) { derived_cast().getHessianProd(d, Hd); }
};

template<class T> struct IsOracle : std::bool_constant< std::is_base_of<OracleBase<T>, T>::value > {};

template<class T> concept Oracle_concept            = IsOracle<T>::value;
template<class T> concept FirstOrderOracle_concept  = IsOracle<T>::value and T::hasGradient;
template<class T> concept SecondOrderOracle_concept = IsOracle<T>::value and T::hasGradient and T::hasHessianProd;

} // namespace LNOT

#endif // LNOT_ORACLE_BASE_HPP

