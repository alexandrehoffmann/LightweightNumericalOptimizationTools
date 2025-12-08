#ifndef LNOT_CRTP_BASE_HPP
#define LNOT_CRTP_BASE_HPP

namespace LNOT
{
	
template<class Derived>
class CRTPBase
{
public:
	      Derived& derived()       { return static_cast<      Derived&>(*this); }
	const Derived& derived() const { return static_cast<const Derived&>(*this); }
};
	
} // namespace FSLinalg

#endif // LNOT_CRTP_BASE_HPP
