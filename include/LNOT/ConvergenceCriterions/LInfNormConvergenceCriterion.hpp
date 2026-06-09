#ifndef LNOT_LInf_NORM_CONVERGENCE_CRITERION_HPP
#define LNOT_LInf_NORM_CONVERGENCE_CRITERION_HPP

#include <LNOT/ConvergenceCriterions/CConvergenceCriterionFor.hpp>

#include <LNOT/BasicLinalg/BasicLinalg.hpp>

namespace LNOT
{

class LInfNormConvergenceCriterion
{
public:
	template<typename Scalar> 
	constexpr inline Scalar getAbsTol(const Scalar& tol) const { return tol; }
	
	template<typename Scalar> 
	constexpr inline Scalar getRelTol(const Scalar& tol, const Scalar& residual) const { return tol*residual; }
	
	template<typename Scalar, typename Size>
	constexpr inline Scalar getResidual(const Scalar* gradient, const Size size) const { return BasicLinalg::lInfNorm(gradient, size); }
};

} // namespace LNOT

#endif // LNOT_LInf_NORM_CONVERGENCE_CRITERION_HPP
