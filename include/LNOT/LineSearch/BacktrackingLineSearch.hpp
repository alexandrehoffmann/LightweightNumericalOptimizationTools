#ifndef LNOT_BACKTRACKING_LINESEARCH_HPP
#define LNOT_BACKTRACKING_LINESEARCH_HPP

#include <LNOT/LineSearch/LineSearchBase.hpp>

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
	
	BacktrackingLineSearch(const Size maxIt = 200000, const Scalar tau = 0.5, const Scalar armijoConditionConst = 1.0e-4) : Base(maxIt), m_tau(tau), m_armijoConditionConst(armijoConditionConst) {}

	void clearWorkSpace();

	template<FirstOrderOracle_concept Oracle>
	Scalar solve(const Scalar* x, const Scalar fx, const Scalar* __restrict__ gx, const Scalar* __restrict__ s, Oracle& oracle);
private:	
	Scalar m_tau;
	Scalar m_armijoConditionConst;
	
	Scalar* m_newX = nullptr;
};
	
} // namespace LNOT

#include <LNOT/LineSearch/BacktrackingLineSearch_impl.hpp>

#endif // LNOT_BACKTRACKING_LINESEARCH_HPP
