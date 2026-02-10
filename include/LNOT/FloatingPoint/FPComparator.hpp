#ifndef LNOT_FP_COMPARATOR_HPP
#define LNOT_FP_COMPARATOR_HPP

#include <cmath>
#include <limits>
#include <algorithm>

namespace LNOT
{

template<class Scalar>
class FPComparator
{
public:
	static inline const Scalar defaultEps = std::numeric_limits<Scalar>::epsilon();

	FPComparator(const Scalar& relativeEpsilon = defaultEps, const Scalar& absoluteEpsilon = defaultEps) : m_relativeEpsilon(relativeEpsilon), m_absoluteEpsilon(absoluteEpsilon) {}
	
	constexpr bool isApproxEq       (const Scalar& a, const Scalar& b) const { return std::abs(a - b) < std::max(std::abs(a), std::abs(b))*m_relativeEpsilon or std::abs(a - b) < m_absoluteEpsilon; }
	constexpr bool isDefGreaterThan (const Scalar& a, const Scalar& b) const { return         (a - b) > std::max(std::abs(a), std::abs(b))*m_relativeEpsilon or         (a - b) > m_absoluteEpsilon; }
	constexpr bool isDefLessThan    (const Scalar& a, const Scalar& b) const { return         (b - a) > std::max(std::abs(a), std::abs(b))*m_relativeEpsilon or         (b - a) > m_absoluteEpsilon; }
	
	constexpr bool isApproxZero  (const Scalar& a) const { return std::abs(a) < m_absoluteEpsilon; }
	constexpr bool isDefPositive (const Scalar& a) const { return          a  > m_absoluteEpsilon; }
	constexpr bool isDefNegative (const Scalar& a) const { return         -a  > m_absoluteEpsilon; }
	
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
