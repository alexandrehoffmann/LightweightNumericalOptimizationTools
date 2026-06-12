#ifndef LNOT_L1_NORM_HPP
#define LNOT_L1_NORM_HPP

#include <LNOT/ConvergenceCriterions/ConvergenceCriterionBase.hpp>

#include <LNOT/BasicLinalg/BasicLinalg.hpp>

namespace LNOT
{

template<typename AScalar> class L1Norm;

template<typename AScalar> 
struct ConvergenceCriterionTraits< L1Norm<AScalar> >
{
	using Size   = unsigned int;
	using Scalar = AScalar;
};

template<typename Scalar_> 
class L1Norm : public ConvergenceCriterionBase< L1Norm<Scalar_> >
{
	using Self = L1Norm<Scalar_>;
public:
	LNOT_DEFINE_CONVERGENCE_CRITERION

	template<typename ASize>
	constexpr void allocate(const ASize /* size */) requires(isSize<ASize>) {}

	constexpr Scalar getTol(const Scalar& tol) const { return tol; }
	
	template<typename ASize>
	constexpr Scalar getResidualImpl(const Scalar* gradient, const ASize size) const requires(isSize<ASize>) { return BasicLinalg::l1Norm(gradient, size); }
};

} // namespace LNOT

#endif // LNOT_L1_NORM_HPP
