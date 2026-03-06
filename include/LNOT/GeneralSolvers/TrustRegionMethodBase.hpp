#ifndef LNOT_TRUST_REGION_METHOD_BASE_HPP
#define LNOT_TRUST_REGION_METHOD_BASE_HPP

#include <LNOT/FloatingPoint/FPComparator.hpp>

namespace LNOT
{

template<typename Scalar, typename Size>
class TrustRegionMethodBase
{
public:
	using Comparator = FPComparator<Scalar>;

	Scalar getEtaVerySuccessful () const { return m_etaVerySuccessful; }
	Scalar getEtaSuccessful     () const { return m_etaSuccessful;     }
	Scalar getEtaAccept         () const { return m_etaAccept;         }
	
	void setEtaVerySuccessful (const Scalar& etaVerySuccessful) { m_etaVerySuccessful = etaVerySuccessful; }
	void setEtaSuccessful     (const Scalar& etaSuccessful)     { m_etaSuccessful     = etaSuccessful;     }
	void setEtaAccept         (const Scalar& etaAccept)         { m_etaAccept         = etaAccept;         }

	Scalar getGammaIncrease() const { return m_gammaIncrease; }
	Scalar getGammaDecrease() const { return m_gammaDecrease; }
	
	void setGammaIncrease(const Scalar& gammaIncrease) { m_gammaIncrease = gammaIncrease; }
	void setGammaDecrease(const Scalar& gammaDecrease) { m_gammaDecrease = gammaDecrease; }
	
	bool isStepSuccessful     (const Scalar& ared, const Scalar& pred, const Comparator& cmp) const { return cmp.isDefGreaterThan(ared, m_etaAccept) and cmp.isDefGreaterThan(ared, pred*m_etaSuccessful);     }
	bool isStepVerySuccessful (const Scalar& ared, const Scalar& pred, const Comparator& cmp) const { return cmp.isDefGreaterThan(ared, m_etaAccept) and cmp.isDefGreaterThan(ared, pred*m_etaVerySuccessful); }
	bool isStepAccepted       (const Scalar& ared, const Scalar& pred, const Comparator& cmp) const { return cmp.isDefGreaterThan(ared, pred*m_etaAccept); }
protected:
	Scalar m_etaVerySuccessful = Scalar(0.75);
	Scalar m_etaSuccessful     = Scalar(0.25);
	Scalar m_etaAccept         = Scalar{};
	
	Scalar m_gammaIncrease = Scalar(2.00);
	Scalar m_gammaDecrease = Scalar(0.25);
	
	Scalar m_deltaMax = Scalar(100.);
};
	
} // namespace LNOT

#define LNOT_DEFINE_TRUST_REGION_SOLVER \
	using TRMBase = TrustRegionMethodBase<Scalar, Size>; \

#define LNOT_TRUST_REGION_SOLVER_ATTRIBUTE \
	using TRMBase::m_etaVerySuccessful; \
	using TRMBase::m_etaSuccessful; \
	using TRMBase::m_etaAccept; \
	using TRMBase::m_gammaIncrease; \
	using TRMBase::m_gammaDecrease; \
	using TRMBase::m_deltaMax; \


#endif // LNOT_TRUST_REGION_METHOD_BASE_HPP
