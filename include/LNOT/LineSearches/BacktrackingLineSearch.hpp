#ifndef LNOT_BACKTRACKING_LINESEARCH_HPP
#define LNOT_BACKTRACKING_LINESEARCH_HPP

#include <LNOT/LineSearches/LineSearchBase.hpp>

namespace LNOT
{
	
template<typename T> class BacktrackingLineSearch;

template<typename T> struct LineSearchTraits< BacktrackingLineSearch<T> > { using Scalar = T; using Size = unsigned int; };

template<typename T> 
class BacktrackingLineSearch : public LineSearchBase< BacktrackingLineSearch<T> >
{
	using Base = LineSearchBase< BacktrackingLineSearch<T> >;
public:
	using Scalar = typename Base::Scalar;
	using Size   = typename Base::Size;
	using Info   = typename Base::Info;
	
	BacktrackingLineSearch(const Size maxIt = 200000, const Scalar tau = Scalar(0.5), const Scalar armijoConditionConst = Scalar(1.0e-4)) : Base(maxIt), m_tau(tau), m_armijoConditionConst(armijoConditionConst) {}

	void clearWorkSpace();

	template<FirstOrderOracle_concept Oracle>
	Scalar solveImpl(const Scalar* x, const Scalar fx, const Scalar* gx, const Scalar* s, Oracle& oracle);
	
	Scalar getTau                   () const { return m_tau;                  }
	Scalar getmArmijoConditionConst () const { return m_armijoConditionConst; }
	
	void setTau                   (const Scalar& tau)                  { m_tau = tau; }
	void setmArmijoConditionConst (const Scalar& armijoConditionConst) { m_armijoConditionConst = armijoConditionConst; }
	
private:	
	Scalar m_tau;
	Scalar m_armijoConditionConst;
	
	Scalar* m_newX = nullptr;
};

extern template class BacktrackingLineSearch<float>;
extern template class BacktrackingLineSearch<double>;
	
} // namespace LNOT

#include <LNOT/LineSearches/BacktrackingLineSearch_impl.hpp>

#endif // LNOT_BACKTRACKING_LINESEARCH_HPP
