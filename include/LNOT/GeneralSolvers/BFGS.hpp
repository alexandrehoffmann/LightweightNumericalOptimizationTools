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
	
	BFGS(LineSearch& lineSearch, const Size maxIt = 200000, const Scalar tol = std::numeric_limits<Scalar>::epsilon()) : Base(maxIt, tol), m_lineSearch(lineSearch) {}
	
	void clearWorkSpaceImpl();
	
	template<FirstOrderOracle_concept Oracle, bool solveInPlace> 
	void solveImpl(Oracle& oracle, std::bool_constant<solveInPlace> bc, Scalar* x);
private:
	Scalar* m_gk    = nullptr;
	Scalar* m_gkp1  = nullptr;
	Scalar* m_sk    = nullptr;
	Scalar* m_yk    = nullptr;
	Scalar* m_uk    = nullptr;
	Scalar* m_invBk = nullptr;
	
	LineSearch& m_lineSearch;
};

template<LineSearch_concept LineSearch>
BFGS<LineSearch> makeBFGS(LineSearch& lineSearch, const typename LineSearch::Size maxIt = 200000, const typename LineSearch::Scalar tol = std::numeric_limits<typename LineSearch::Scalar>::epsilon()) { return BFGS<LineSearch>(lineSearch, maxIt, tol); }

} // namespace LNOT

#include <LNOT/GeneralSolvers/BFGS_impl.hpp>

#endif // LNOT_BFGS_HPP
