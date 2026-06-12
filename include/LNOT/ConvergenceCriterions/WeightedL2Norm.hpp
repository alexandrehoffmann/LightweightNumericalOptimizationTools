#ifndef LNOT_WEIGHTED_L2_NORM_HPP
#define LNOT_WEIGHTED_L2_NORM_HPP

#include <LNOT/ConvergenceCriterions/ConvergenceCriterionBase.hpp>

#include <LNOT/BasicLinalg/BasicLinalg.hpp>

#include <vector>
#include <span>
#include <cassert>

namespace LNOT
{

template<typename AScalar> class WeightedL2Norm;

template<typename AScalar> 
struct ConvergenceCriterionTraits< WeightedL2Norm<AScalar> >
{
	using Size   = unsigned int;
	using Scalar = AScalar;
};

template<typename Scalar_> 
class WeightedL2Norm : public ConvergenceCriterionBase< WeightedL2Norm<Scalar_> >
{
	using Self = WeightedL2Norm<Scalar_>;
public:
	LNOT_DEFINE_CONVERGENCE_CRITERION
	
	template<typename ASize>
	void allocate(const ASize size) requires(isSize<ASize>) { m_weights.resize(size); }

	constexpr std::span<const Scalar> getWeights() const { return std::span(m_weights.get(), m_size); }
	
	constexpr std::span<Scalar> getWeights() { return std::span(m_weights.get(), m_size); }

	constexpr Scalar getTol(const Scalar& tol) const { return tol*tol; }
	
	template<typename ASize>
	constexpr Scalar getResidualImpl(const Scalar* gradient, const ASize size) const requires(isSize<ASize>) { assert(size == m_size); return BasicLinalg::weightedSquaredNorm(gradient, m_weights.get(), size); }
protected:
	std::vector<Scalar> m_weights;
};

} // namespace LNOT

#endif // LNOT_WEIGHTED_L2_NORM_HPP
