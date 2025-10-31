#ifndef LNOT_LBFGS_HPP
#define LNOT_LBFGS_HPP

#include <LNOT/GeneralSolvers/FirstOrderSolverBase.hpp>
#include <LNOT/GeneralSolvers/LimitedMemorySolverBase.hpp>
#include <LNOT/LineSearches/LineSearchBase.hpp>

namespace LNOT
{

template<typename LineSearch> class LBFGS;

template<typename LineSearch> 
struct FirstOrderSolverTraits< LBFGS<LineSearch> >
{
	static_assert(LineSearch_concept<LineSearch>);
	
	using Scalar = typename LineSearch::Scalar;
	using Size   = typename LineSearch::Size;
};

template<typename LineSearch>
class LBFGS 
	: public FirstOrderSolverBase< LBFGS<LineSearch> >
	, public LimitedMemorySolverBase<typename LineSearch::Size>
{
	using Base   = FirstOrderSolverBase< LBFGS<LineSearch> >;
	using LMBase = LimitedMemorySolverBase<typename LineSearch::Size>;
public:
	using Scalar = typename Base::Scalar;
	using Size   = typename Base::Size;
	using Info   = typename Base::Info;
	
	void clearWorkSpaceImpl();
	
	template<FirstOrderOracle_concept Oracle, bool solveInPlace> 
	void solveImpl(Oracle& oracle, std::bool_constant<solveInPlace> bc, Scalar* x);
	
	const LineSearch& getLinesearch() const { return m_lineSearch; }
	      LineSearch& getLinesearch()       { return m_lineSearch; }
private:
	Scalar* m_gk   = nullptr;
	Scalar* m_gkp1 = nullptr;
	Scalar* m_S    = nullptr;
	Scalar* m_Y    = nullptr;
	Scalar* m_dk   = nullptr;
	
	LineSearch m_lineSearch;
};

} // namespace LNOT

#include <LNOT/GeneralSolvers/LBFGS_impl.hpp>

#endif // LNOT_LBFGS_HPP
