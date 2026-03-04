#ifndef LNOT_FP_COMPARATOR_HPP
#define LNOT_FP_COMPARATOR_HPP

#include <LNOT/FloatingPoint/NumTraits.hpp>
#include <LNOT/misc/AdlMath.hpp>

#include <cmath>
#include <algorithm>

namespace LNOT
{

template<std::floating_point Scalar>
class FPComparator
{
public:
	static inline constexpr Scalar defaultEps = NumTraits<Scalar>::epsilon;

	FPComparator(const Scalar& relativeEpsilon = defaultEps, const Scalar& absoluteEpsilon = defaultEps) : m_relativeEpsilon(relativeEpsilon), m_absoluteEpsilon(absoluteEpsilon) {}
	
	constexpr bool isApproxEq(const Scalar& a, const Scalar& b) const { return AdlMath::abs(a - b) < std::max(AdlMath::abs(a), AdlMath::abs(b))*m_relativeEpsilon or AdlMath::abs(a - b) < m_absoluteEpsilon; }
	
	constexpr bool isDefGreaterThan(const Scalar& a, const Scalar& b) const { return (a - b) > std::max(AdlMath::abs(a), AdlMath::abs(b))*m_relativeEpsilon or (a - b) > m_absoluteEpsilon; }
	
	constexpr bool isDefLessThan(const Scalar& a, const Scalar& b) const { return (b - a) > std::max(AdlMath::abs(a), AdlMath::abs(b))*m_relativeEpsilon or (b - a) > m_absoluteEpsilon; }
	
	constexpr bool isApproxZero(const Scalar& a) const { return AdlMath::abs(a) < m_absoluteEpsilon; }
	
	constexpr bool isDefPositive(const Scalar& a) const { return a  > m_absoluteEpsilon; }
	
	constexpr bool isDefNegative(const Scalar& a) const { return -a  > m_absoluteEpsilon; }
	
	constexpr Scalar getRelativeEpsilon() const { return m_relativeEpsilon; }
	
	constexpr Scalar getAbsoluteEpsilon() const { return m_absoluteEpsilon; }
	
	void setRelativeEpsilon(const Scalar& epsilon) const { m_relativeEpsilon = epsilon; }
	
	void setAbsoluteEpsilon(const Scalar& epsilon) const { m_absoluteEpsilon = epsilon; }
private:
	Scalar m_relativeEpsilon;
	Scalar m_absoluteEpsilon;
};

} // namespace LNOT

#endif // LNOT_FP_COMPARATOR_HPP
