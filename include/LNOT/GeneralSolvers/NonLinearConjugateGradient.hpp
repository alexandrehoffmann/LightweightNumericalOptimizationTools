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
	static_assert(LineSearch_concept<LineSearch>);
	
	using Scalar = typename LineSearch::Scalar;
	using Size   = typename LineSearch::Size;
};

template<typename LineSearch, NLCGUpdateStrategy UpdateStrategy = NLCGUpdateStrategy::HESTENES_STIEFEL> 
class NonLinearConjugateGradient : public FirstOrderSolverBase< NonLinearConjugateGradient<LineSearch, UpdateStrategy> >
{
	using Base = FirstOrderSolverBase< NonLinearConjugateGradient<LineSearch, UpdateStrategy> >;
public:
	using Scalar = typename Base::Scalar;
	using Size   = typename Base::Size;
	using Info   = typename Base::Info;
	
	NonLinearConjugateGradient(LineSearch& lineSearch, const Size maxIt = 200000, const Scalar tol = std::numeric_limits<Scalar>::epsilon()) : Base(maxIt, tol), m_lineSearch(lineSearch) {}
	
	void clearWorkSpaceImpl();
	
	template<FirstOrderOracle_concept Oracle, bool solveInPlace> 
	void solveImpl(Oracle& oracle, std::bool_constant<solveInPlace> bc, Scalar* x);
	
	Scalar getBeta(const Size size);
private:
	Scalar* m_dk   = nullptr;
	Scalar* m_yk   = nullptr;
	Scalar* m_gk   = nullptr;
	Scalar* m_gkp1 = nullptr;
	
	LineSearch& m_lineSearch;
};


template<LineSearch_concept LineSearch, NLCGUpdateStrategy UpdateStrategy>
NonLinearConjugateGradient<LineSearch, UpdateStrategy> makeNLCG(LineSearch& lineSearch, const NLCGUpdateStrategy_constant<UpdateStrategy>, const typename LineSearch::Size maxIt = 200000, const typename LineSearch::Scalar tol = std::numeric_limits<typename LineSearch::Scalar>::epsilon()) { return NonLinearConjugateGradient<LineSearch, UpdateStrategy>(lineSearch, maxIt, tol); }

template<LineSearch_concept LineSearch>
NonLinearConjugateGradient<LineSearch> makeNLCG(LineSearch& lineSearch, const typename LineSearch::Size maxIt = 200000, const typename LineSearch::Scalar tol = std::numeric_limits<typename LineSearch::Scalar>::epsilon()) { return NonLinearConjugateGradient<LineSearch>(lineSearch, maxIt, tol); }


} // namespace LNOT

#include <LNOT/GeneralSolvers/NonLinearConjugateGradient_impl.hpp>

#endif // LNOT_NON_LINEAR_CONJUGATE_GRADIENT_HPP
