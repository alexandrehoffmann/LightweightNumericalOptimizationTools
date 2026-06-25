#ifndef LNOT_CONVERGENCE_CRITERION_BASE_HPP
#define LNOT_CONVERGENCE_CRITERION_BASE_HPP

#include <LNOT/CRTPBase.hpp>

namespace LNOT
{

template<class Derived> struct ConvergenceCriterionTraits;

template<class Derived>
class ConvergenceCriterionBase : public CRTPBase<Derived>
{
	using DerivedTraits = ConvergenceCriterionTraits<Derived>;
public:
	using Size   = typename DerivedTraits::Size;   ///<  @brief Type representing the number of dimensions (variables) of the function.
	using Scalar = typename DerivedTraits::Scalar; ///<  @brief Scalar type used in evaluations.

	template<typename ASize> static constexpr bool isSize = std::same_as<Size, BIC::Mutable<ASize>>;
	
	template<typename ASize>
	constexpr void init(const Scalar* grad, const ASize size) requires(isSize<ASize>) { this->derived().allocate(size); m_initialResidual = getResidual(grad, size); }
	
	constexpr const Scalar& getInitialResidual() const { return m_initialResidual; }
	
	template<typename ASize>
	constexpr Scalar getResidual(const Scalar* grad, const ASize size) const requires(isSize<ASize>) { return this->derived().getResidualImpl(grad, size); } 
	
	constexpr std::convertible_to<Scalar> auto getAbsTol(const Scalar& tol) const { return this->derived().getTol(tol); }
	
	constexpr Scalar getRelTol(const Scalar& tol) const { return m_initialResidual*this->derived().getTol(tol); }
private:
	Scalar m_initialResidual = {};
};

#define LNOT_DEFINE_CONVERGENCE_CRITERION \
	using Base   = LNOT::ConvergenceCriterionBase<Self>; \
	using Size   = typename Base::Size; \
	using Scalar = typename Base::Scalar; \
	\
	template<typename ASize> static constexpr bool isSize = Base::template isSize<ASize>; \

template<class T> struct IsConvergenceCriterion : BIC::Fixed<bool, std::is_base_of<ConvergenceCriterionBase<T>, T>::value > {};

template<class T> concept CConvergenceCriterion = IsConvergenceCriterion<T>::value;

} // namespace LNOT

#endif // LNOT_CONVERGENCE_CRITERION_BASE_HPP
