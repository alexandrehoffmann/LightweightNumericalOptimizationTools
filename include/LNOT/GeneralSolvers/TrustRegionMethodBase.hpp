#ifndef LNOT_TRUST_REGION_METHOD_BASE_HPP
#define LNOT_TRUST_REGION_METHOD_BASE_HPP

namespace LNOT
{

template<typename Scalar, typename Size>
class TrustRegionMethodBase
{
public:
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
protected:
	Scalar m_etaVerySuccessful = Scalar(0.75);
	Scalar m_etaSuccessful     = Scalar(0.25);
	Scalar m_etaAccept         = Scalar(0.10);
	
	Scalar m_gammaIncrease = Scalar(2.00);
	Scalar m_gammaDecrease = Scalar(0.25);
};
	
}

#endif // LNOT_TRUST_REGION_METHOD_BASE_HPP
