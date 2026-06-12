#ifndef LNOT_WEIGHTED_L2_NORM_IMPL_HPP
#define LNOT_WEIGHTED_L2_NORM_IMPL_HPP

namespace LNOT
{

template<typename Scalar_> template<typename ASize> 
void WeightedL2Norm<Scalar_>::allocate(const ASize size) requires(isSize<ASize>)
{
	m_size = size;
	if (size > m_capacity) 
	{ 
		m_capacity = size;		
		m_weights = std::make_unique<Scalar[]>(m_capacity);
	}
}

} // namespace LNOT

#endif // LNOT_WEIGHTED_L2_NORM_IMPL_HPP
