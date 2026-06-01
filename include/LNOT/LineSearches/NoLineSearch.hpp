#ifndef LNOT_NO_LINE_SEARCH_HPP
#define LNOT_NO_LINE_SEARCH_HPP

#include <LNOT/LineSearches/LineSearchBase.hpp>

namespace LNOT
{

template<typename T> class NoLineSearch;


template<typename T> struct LineSearchTraits< NoLineSearch<T> > { using Scalar = T; using Size = unsigned int; };

template<typename T>
class NoLineSearch : public LineSearchBase< NoLineSearch<T> >
{
	using Self = NoLineSearch<T>;
public:
	LNOT_DEFINE_LINESEARCH
	
	NoLineSearch(const Scalar alpha = Scalar(1)) : Base(1), m_alpha(alpha) {}
	
	template<CFirstOrderOracle Oracle>
	Scalar solveImpl(const Scalar* /* x */, const Scalar& /* fx */, const Scalar* /* gx */, const Scalar* /* s */, Oracle& /* oracle */) { return m_alpha; }
	
	Scalar getAlpha() const { return m_alpha; }
	
	void setAlpha(const Scalar alpha) { m_alpha = alpha; }
protected:
	LNOT_LINESEARCH_ATTRIBUTE
private:
	Scalar m_alpha;
};

extern template class NoLineSearch<float>;
extern template class NoLineSearch<double>;
extern template class NoLineSearch<long double>;

} // namespace LNOT

#endif // LNOT_NO_LINE_SEARCH_HPP
