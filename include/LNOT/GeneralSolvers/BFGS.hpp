#ifndef LNOT_BFGS_HPP
#define LNOT_BFGS_HPP

#include <LNOT/GeneralSolvers/FirstOrderSolverBase.hpp>
#include <LNOT/LineSearches/LineSearchBase.hpp>

#include <memory>

namespace LNOT
{

template<typename LineSearch> class BFGS;

template<typename LineSearch> 
struct FirstOrderSolverTraits< BFGS<LineSearch> >
{
	static_assert(CLineSearch<LineSearch>);
	
	using Scalar = typename LineSearch::Scalar;
	using Size   = typename LineSearch::Size;
};

template<typename LineSearch>
class BFGS : public FirstOrderSolverBase< BFGS<LineSearch> >
{
	using Self = BFGS<LineSearch>;
public:
	LNOT_DEFINE_FIRST_ORDER_SOLVER
	
	template<CFirstOrderOracle Oracle, typename ABool> 
	void solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value);
	
	const LineSearch& getLinesearch() const { return m_lineSearch; }
	      LineSearch& getLinesearch()       { return m_lineSearch; }
protected:
	LNOT_FIRST_ORDER_SOLVER_ATTRIBUTE
private:
	std::unique_ptr<Scalar[]> m_gk;
	std::unique_ptr<Scalar[]> m_gkp1;
	std::unique_ptr<Scalar[]> m_sk;
	std::unique_ptr<Scalar[]> m_yk;
	std::unique_ptr<Scalar[]> m_uk;
	std::unique_ptr<Scalar[]> m_invBk;
	
	LineSearch m_lineSearch;
};

} // namespace LNOT

#include <LNOT/GeneralSolvers/BFGS_impl.hpp>

#endif // LNOT_BFGS_HPP
