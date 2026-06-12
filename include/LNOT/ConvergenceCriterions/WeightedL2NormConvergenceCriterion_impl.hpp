#ifndef LNOT_WEIGHTED_L2_NORM_CONVERGENCE_CRITERION_IMPL_HPP
#define LNOT_WEIGHTED_L2_NORM_CONVERGENCE_CRITERION_IMPL_HPP

namespace LNOT
{

template<typename Scalar, typename Size>
WeightedL2NormConvergenceCriterion::Scope<Scalar,Size>::Scope(const Size size, const Scalar& value) 
	: m_previousWeights(std::move(s_activeWeights))
	, m_previousSize(std::move(s_activeSize))
{ 
	s_activeWeights = std::make_unique<Scalar[]>(size);
	s_activeSize    = size;
	
	std::fill(s_activeWeights.get(), s_activeWeights.get() + size, value);
}
		
template<typename Scalar, typename Size>
WeightedL2NormConvergenceCriterion::Scope<Scalar,Size>::Scope(const Scalar* values, const Size size)
	: m_previousWeights(std::move(s_activeWeights))
	, m_previousSize(std::move(s_activeSize))
{ 
	s_activeWeights = std::make_unique<Scalar[]>(size);
	s_activeSize    = size;
	
	std::copy(values, values + size, s_activeWeights.get());
}
		
template<typename Scalar, typename Size>
WeightedL2NormConvergenceCriterion::Scope<Scalar,Size>::~Scope()
{
	s_activeWeights = std::move(m_previousWeights);
	s_activeSize    = std::move(m_previousSize);
}

} // namespace LNOT

#endif // LNOT_WEIGHTED_L2_NORM_CONVERGENCE_CRITERION_IMPL_HPP
