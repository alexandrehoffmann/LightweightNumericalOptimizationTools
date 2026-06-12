#ifndef LNOT_L2_NORM_HPP
#define LNOT_L2_NORM_HPP

#include <LNOT/ConvergenceCriterions/ConvergenceCriterionBase.hpp>

#include <LNOT/BasicLinalg/BasicLinalg.hpp>

namespace LNOT
{

template<typename AScalar> class L2Norm;

template<typename AScalar> 
struct ConvergenceCriterionTraits< L2Norm<AScalar> >
{
	using Size   = unsigned int;
	using Scalar = AScalar;
};

template<typename Scalar_> 
class L2Norm : public ConvergenceCriterionBase< L2Norm<Scalar_> >
{
	using Self = L2Norm<Scalar_>;
public:
	LNOT_DEFINE_CONVERGENCE_CRITERION
	
	template<typename ASize>
	constexpr void allocate(const ASize /* size */) requires(isSize<ASize>) {}

	constexpr Scalar getTol(const Scalar& tol) const { return tol*tol; }
	
	template<typename ASize>
	constexpr Scalar getResidualImpl(const Scalar* gradient, const ASize size) const requires(isSize<ASize>) { return BasicLinalg::squaredNorm(gradient, size); }
};

} // namespace LNOT

#endif // LNOT_L2_NORM_HPP
