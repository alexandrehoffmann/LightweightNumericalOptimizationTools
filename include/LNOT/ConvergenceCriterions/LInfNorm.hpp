#ifndef LNOT_LINF_NORM_HPP
#define LNOT_LINF_NORM_HPP

#include <LNOT/ConvergenceCriterions/ConvergenceCriterionBase.hpp>

#include <LNOT/BasicLinalg/BasicLinalg.hpp>

namespace LNOT
{

template<typename AScalar> class LInfNorm;

template<typename AScalar> 
struct ConvergenceCriterionTraits< LInfNorm<AScalar> >
{
	using Size   = unsigned int;
	using Scalar = AScalar;
};

template<typename Scalar_> 
class LInfNorm : public ConvergenceCriterionBase< LInfNorm<Scalar_> >
{
	using Self = LInfNorm<Scalar_>;
public:
	LNOT_DEFINE_CONVERGENCE_CRITERION
	
	template<typename ASize>
	constexpr void allocate(const ASize /* size */) requires(isSize<ASize>) {}

	constexpr const Scalar& getTol(const Scalar& tol) const { return tol; }
	
	template<typename ASize>
	constexpr Scalar getResidualImpl(const Scalar* gradient, const ASize size) const requires(isSize<ASize>) { return BasicLinalg::lInfNorm(gradient, size); }
};

} // namespace LNOT

#endif // LNOT_LINF_NORM_HPP
