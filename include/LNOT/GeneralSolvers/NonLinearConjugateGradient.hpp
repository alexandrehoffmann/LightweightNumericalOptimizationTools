#ifndef LNOT_NON_LINEAR_CONJUGATE_GRADIENT_HPP
#define LNOT_NON_LINEAR_CONJUGATE_GRADIENT_HPP

#include <LNOT/GeneralSolvers/FirstOrderSolverBase.hpp>
#include <LNOT/GeneralSolvers/NLCGUpdateStrategy.hpp>
#include <LNOT/LineSearches/LineSearchBase.hpp>

namespace LNOT
{

template<typename LineSearch, NLCGUpdateStrategy UpdateStrategy> class NonLinearConjugateGradient;

template<typename LineSearch, NLCGUpdateStrategy UpdateStrategy> 
struct FirstOrderSolverTraits< NonLinearConjugateGradient<LineSearch, UpdateStrategy> >
{
	static_assert(CLineSearch<LineSearch>);
	
	using Scalar = typename LineSearch::Scalar;
	using Size   = typename LineSearch::Size;
};

template<typename LineSearch, NLCGUpdateStrategy UpdateStrategy = NLCGUpdateStrategy::HESTENES_STIEFEL> 
class NonLinearConjugateGradient : public FirstOrderSolverBase< NonLinearConjugateGradient<LineSearch, UpdateStrategy> >
{
	using Self = NonLinearConjugateGradient<LineSearch, UpdateStrategy>;
public:
	LNOT_DEFINE_FIRST_ORDER_SOLVER
	
	void clearWorkSpaceImpl();
	
	template<CFirstOrderOracle Oracle, typename ABool> 
	void solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value);
	
	Scalar getBeta(const Size size);
	
	const LineSearch& getLinesearch() const { return m_lineSearch; }
	      LineSearch& getLinesearch()       { return m_lineSearch; }
protected:
	LNOT_FIRST_ORDER_SOLVER_ATTRIBUTE
private:
	Scalar* m_dk   = nullptr;
	Scalar* m_yk   = nullptr;
	Scalar* m_gk   = nullptr;
	Scalar* m_gkp1 = nullptr;
	
	LineSearch m_lineSearch;
};

} // namespace LNOT

#include <LNOT/GeneralSolvers/NonLinearConjugateGradient_impl.hpp>

#endif // LNOT_NON_LINEAR_CONJUGATE_GRADIENT_HPP
