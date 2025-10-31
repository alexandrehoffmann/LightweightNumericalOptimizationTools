#ifndef LNOT_BFGS_HPP
#define LNOT_BFGS_HPP

#include <LNOT/GeneralSolvers/FirstOrderSolverBase.hpp>
#include <LNOT/LineSearches/LineSearchBase.hpp>

namespace LNOT
{

template<typename LineSearch> class BFGS;

template<typename LineSearch> 
struct FirstOrderSolverTraits< BFGS<LineSearch> >
{
	static_assert(LineSearch_concept<LineSearch>);
	
	using Scalar = typename LineSearch::Scalar;
	using Size   = typename LineSearch::Size;
};

template<typename LineSearch>
class BFGS : public FirstOrderSolverBase< BFGS<LineSearch> >
{
	using Base = FirstOrderSolverBase< BFGS<LineSearch> >;
public:
	using Scalar = typename Base::Scalar;
	using Size   = typename Base::Size;
	using Info   = typename Base::Info;
	
	void clearWorkSpaceImpl();
	
	template<FirstOrderOracle_concept Oracle, bool solveInPlace> 
	void solveImpl(Oracle& oracle, std::bool_constant<solveInPlace>, Scalar* x);
	
	const LineSearch& getLinesearch() const { return m_lineSearch; }
	      LineSearch& getLinesearch()       { return m_lineSearch; }
private:
	Scalar* m_gk    = nullptr;
	Scalar* m_gkp1  = nullptr;
	Scalar* m_sk    = nullptr;
	Scalar* m_yk    = nullptr;
	Scalar* m_uk    = nullptr;
	Scalar* m_invBk = nullptr;
	
	LineSearch m_lineSearch;
};

} // namespace LNOT

#include <LNOT/GeneralSolvers/BFGS_impl.hpp>

#endif // LNOT_BFGS_HPP
