#ifndef LNOT_LBFGS_HPP
#define LNOT_LBFGS_HPP

#include <LNOT/GeneralSolvers/FirstOrderSolverBase.hpp>

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
class LBFGS : public FirstOrderSolverBase< LBFGS<LineSearch> >
{
	using Base = FirstOrderSolverBase< LBFGS<LineSearch> >;
public:
	using Scalar = typename Base::Scalar;
	using Size   = typename Base::Size;
	using Info   = typename Base::Info;
	
	LBFGS(LineSearch& lineSearch, const Size memory=10, const Size maxIt = 200000, const Scalar tol = std::numeric_limits<Scalar>::epsilon()) : Base(maxIt, tol), m_memory(memory), m_lineSearch(lineSearch) {}
	
	void clearWorkSpace();
	
	template<FirstOrderOracle_concept Oracle, bool solveInPlace> 
	void solve_impl(Oracle& oracle, std::bool_constant<solveInPlace> bc, Scalar* x);
private:
	Scalar* m_gk   = nullptr;
	Scalar* m_gkp1 = nullptr;
	Scalar* m_S    = nullptr;
	Scalar* m_Y    = nullptr;
	Scalar* m_dk   = nullptr;
	
	Size m_memory;
	
	LineSearch& m_lineSearch;
};

template<LineSearch_concept LineSearch>
LBFGS<LineSearch> makeLBFGS(LineSearch& lineSearch, const typename LineSearch::Size memory, const typename LineSearch::Size maxIt = 200000, const typename LineSearch::Scalar tol = std::numeric_limits<typename LineSearch::Scalar>::epsilon()) { return LBFGS<LineSearch>(lineSearch, memory, maxIt, tol); }

} // namespace LNOT

#include <LNOT/GeneralSolvers/LBFGS_impl.hpp>

#endif // LNOT_LBFGS_HPP
