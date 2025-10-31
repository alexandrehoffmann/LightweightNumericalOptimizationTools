#ifndef LNOT_BISECTION_LINE_SEARCH_HPP
#define LNOT_BISECTION_LINE_SEARCH_HPP

#include <LNOT/LineSearches/LineSearchBase.hpp>

namespace LNOT
{

template<typename T> class BisectionLineSearch;

template<typename T> struct LineSearchTraits< BisectionLineSearch<T> > { using Scalar = T; using Size = unsigned int; };

template<typename T>
class BisectionLineSearch : public LineSearchBase< BisectionLineSearch<T> >
{
	using Base = LineSearchBase< BisectionLineSearch<T> >;
public:
	using Scalar = typename Base::Scalar;
	using Size   = typename Base::Size;
	using Info   = typename Base::Info;
	
	BisectionLineSearch(const Size maxIt = 200000, const Scalar firstWolfConditionConst = Scalar(1.0e-4), const Scalar secondWolfConditionConst = Scalar(0.9));
	
	void clearWorkSpace();
	
	template<FirstOrderOracle_concept Oracle>
	Scalar solveImpl(const Scalar* x, const Scalar& fx, const Scalar* gx, const Scalar* s, Oracle& oracle);
	
	Scalar getFirstWolfConditionConst  () const { return m_firstWolfConditionConst;  }
	Scalar getSecondWolfConditionConst () const { return m_secondWolfConditionConst; }
	
	void setFirstWolfConditionConst  (const Scalar& firstWolfConditionConst)  { m_firstWolfConditionConst  = firstWolfConditionConst;  }
	void setSecondWolfConditionConst (const Scalar& secondWolfConditionConst) { m_secondWolfConditionConst = secondWolfConditionConst; }
	
private:
	Scalar m_firstWolfConditionConst;
	Scalar m_secondWolfConditionConst;

	Scalar* m_newX    = nullptr;
	Scalar* m_newGrad = nullptr;
};

extern template class BisectionLineSearch<float>;
extern template class BisectionLineSearch<double>;

}

#include <LNOT/LineSearches/BisectionLineSearch_impl.hpp>

#endif // LNOT_BISECTION_LINE_SEARCH_HPP
