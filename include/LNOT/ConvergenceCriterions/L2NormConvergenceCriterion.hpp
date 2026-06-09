#ifndef LNOT_L2_NORM_CONVERGENCE_CRITERION_HPP
#define LNOT_L2_NORM_CONVERGENCE_CRITERION_HPP

#include <LNOT/ConvergenceCriterions/CConvergenceCriterionFor.hpp>

#include <LNOT/BasicLinalg/BasicLinalg.hpp>

namespace LNOT
{

class L2NormConvergenceCriterion
{
public:
	template<typename Scalar> 
	constexpr inline Scalar getAbsTol(const Scalar& tol) const { return tol*tol; }
	
	template<typename Scalar> 
	constexpr inline Scalar getRelTol(const Scalar& tol, const Scalar& residual) const { return tol*tol*residual; }
	
	template<typename Scalar, typename Size>
	constexpr inline Scalar getResidual(const Scalar* gradient, const Size size) const { return BasicLinalg::squaredNorm(gradient, size); }
};

} // namespace LNOT

#endif // LNOT_L2_NORM_CONVERGENCE_CRITERION_HPP
