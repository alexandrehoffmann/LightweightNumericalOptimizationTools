#ifndef LNOT_NON_LINEAR_CONJUGATE_GRADIENT_IMPL_HPP
#define LNOT_NON_LINEAR_CONJUGATE_GRADIENT_IMPL_HPP

#include <LNOT/GeneralSolvers/NonLinearConjugateGradient.hpp>
#include <LNOT/LineSearches/BisectionLineSearch.hpp>
#include <LNOT/LineSearches/BacktrackingLineSearch.hpp>

namespace LNOT
{
	
//// explicit template instanciations ////

//// for floats
extern template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::HESTENES_STIEFEL>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::FLETCHER_REEVES>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::POLAK_RIBIERE>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::FLETCHER>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::LIU_STOREY>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::DAI_YUAN>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::HAGER_ZHANG>;

extern template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::HESTENES_STIEFEL>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::FLETCHER_REEVES>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::POLAK_RIBIERE>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::FLETCHER>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::LIU_STOREY>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::DAI_YUAN>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::HAGER_ZHANG>;

//// for doubles
extern template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::HESTENES_STIEFEL>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::FLETCHER_REEVES>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::POLAK_RIBIERE>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::FLETCHER>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::LIU_STOREY>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::DAI_YUAN>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::HAGER_ZHANG>;

extern template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::HESTENES_STIEFEL>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::FLETCHER_REEVES>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::POLAK_RIBIERE>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::FLETCHER>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::LIU_STOREY>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::DAI_YUAN>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::HAGER_ZHANG>;

//// for long doubles
extern template class NonLinearConjugateGradient<BisectionLineSearch<long double>, NLCGUpdateStrategy::HESTENES_STIEFEL>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<long double>, NLCGUpdateStrategy::FLETCHER_REEVES>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<long double>, NLCGUpdateStrategy::POLAK_RIBIERE>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<long double>, NLCGUpdateStrategy::FLETCHER>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<long double>, NLCGUpdateStrategy::LIU_STOREY>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<long double>, NLCGUpdateStrategy::DAI_YUAN>;
extern template class NonLinearConjugateGradient<BisectionLineSearch<long double>, NLCGUpdateStrategy::HAGER_ZHANG>;

extern template class NonLinearConjugateGradient<BacktrackingLineSearch<long double>, NLCGUpdateStrategy::HESTENES_STIEFEL>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<long double>, NLCGUpdateStrategy::FLETCHER_REEVES>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<long double>, NLCGUpdateStrategy::POLAK_RIBIERE>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<long double>, NLCGUpdateStrategy::FLETCHER>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<long double>, NLCGUpdateStrategy::LIU_STOREY>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<long double>, NLCGUpdateStrategy::DAI_YUAN>;
extern template class NonLinearConjugateGradient<BacktrackingLineSearch<long double>, NLCGUpdateStrategy::HAGER_ZHANG>;

//// method implementations ////

template<typename LineSearch, NLCGUpdateStrategy UpdateStrategy>
void NonLinearConjugateGradient<LineSearch, UpdateStrategy>::clearWorkSpaceImpl()
{
	if (m_dk   != nullptr) { delete[] m_dk;   m_dk   = nullptr; }
	if (m_yk   != nullptr) { delete[] m_yk;   m_yk   = nullptr; }
	if (m_gk   != nullptr) { delete[] m_gk;   m_gk   = nullptr; }
	if (m_gkp1 != nullptr) { delete[] m_gkp1; m_gkp1 = nullptr; }
	Base::m_workCapacity = 0;	
}

template<typename LineSearch, NLCGUpdateStrategy UpdateStrategy>  template<FirstOrderOracle_concept Oracle, bool solveInPlace> 
void NonLinearConjugateGradient<LineSearch, UpdateStrategy>::solveImpl(Oracle& oracle, std::bool_constant<solveInPlace>, Scalar* x)
{
	using Oracle_Size = typename Oracle::Size;
	
	const Oracle_Size size = oracle.getNDims();
	
	if (Base::m_workCapacity < size)
	{
		clearWorkSpaceImpl();
		Base::m_workCapacity = size;

		m_dk   = new Scalar[Base::m_workCapacity];
		m_yk   = new Scalar[Base::m_workCapacity];
		m_gk   = new Scalar[Base::m_workCapacity];
		m_gkp1 = new Scalar[Base::m_workCapacity];
	}
	if constexpr (not solveInPlace) { std::fill(x, x + size, 0); }
	
	Base::m_innerIts.clear();
	
	oracle.setCurrentPoint(x);
	oracle.getGradient(m_gk);
	
	#pragma omp simd
	for (Size i=0; i!=size; ++i) { m_dk[i] = -m_gk[i]; }
	
	Base::m_fx = oracle.getValue();
	Base::m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk, size);
	
	const Scalar tol2 = Base::m_tol*Base::m_tol*std::max(Scalar(1), Base::m_squaredNormGrad);
	
	if (Base::m_out != nullptr) { fmt::print(Base::m_out, "#Non Linear CG method\n#Iteration f(x) residual tol\n"); }
	
	Base::m_info = Info::FAILURE;
	for (Base::m_nIt=0;Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e} {:10.2e}\n", Base::m_nIt, Base::m_fx, Base::m_squaredNormGrad, tol2); }
		if (Base::m_squaredNormGrad < tol2) { Base::m_info = Info::SUCCESS; break; }
		
		const Scalar alpha = m_lineSearch.solve(x, Base::m_fx, m_gk, m_dk, oracle);
		
		if (alpha < std::numeric_limits<Scalar>::epsilon()) { Base::m_info = Info::BREAKDOWN; break; }
		
		Base::m_innerIts.push_back(1);
		
		BasicLinalg::axpy(alpha, m_dk, size, x);
		
		oracle.setCurrentPoint(x);
		oracle.getGradient(m_gkp1);
		
		Base::m_fx = oracle.getValue();
		Base::m_squaredNormGrad = BasicLinalg::squaredNorm(m_gkp1, size);
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_yk[i] = m_gkp1[i] - m_gk[i]; }
		
		const Scalar beta = getBeta(size);
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_dk[i] = -m_gkp1[i] + beta*m_dk[i]; }
		
		std::copy(m_gkp1, m_gkp1 + size, m_gk);
	}
}

template<typename LineSearch, NLCGUpdateStrategy UpdateStrategy> 
auto NonLinearConjugateGradient<LineSearch, UpdateStrategy>::getBeta(const Size size) -> Scalar
{
	if      constexpr (UpdateStrategy == NLCGUpdateStrategy::HESTENES_STIEFEL) { return BasicLinalg::inner(m_gkp1, m_yk, size) /  BasicLinalg::inner(m_dk, m_yk, size); }
	else if constexpr (UpdateStrategy == NLCGUpdateStrategy::FLETCHER_REEVES)  { return BasicLinalg::squaredNorm(m_gkp1, size) /  BasicLinalg::squaredNorm(m_gk, size); }
	else if constexpr (UpdateStrategy == NLCGUpdateStrategy::POLAK_RIBIERE)    { return BasicLinalg::inner(m_gkp1, m_yk, size) /  BasicLinalg::squaredNorm(m_gk, size); }
	else if constexpr (UpdateStrategy == NLCGUpdateStrategy::FLETCHER)         { return BasicLinalg::squaredNorm(m_gkp1, size) / -BasicLinalg::inner(m_dk, m_gk, size); }
	else if constexpr (UpdateStrategy == NLCGUpdateStrategy::LIU_STOREY)       { return BasicLinalg::inner(m_gkp1, m_yk, size) / -BasicLinalg::inner(m_dk, m_gk, size); }
	else if constexpr (UpdateStrategy == NLCGUpdateStrategy::DAI_YUAN)         { return BasicLinalg::squaredNorm(m_gkp1, size) /  BasicLinalg::inner(m_dk, m_yk, size); }
	else                                                                       
	{ 
		const Scalar sqNormY = BasicLinalg::squaredNorm(m_yk, size);
		const Scalar yDotG = BasicLinalg::inner(m_dk, m_gk, size);
		const Scalar beta1 = BasicLinalg::inner(m_yk, m_gkp1, size) / yDotG;
		const Scalar beta2 = -2*(sqNormY / (yDotG*yDotG))*BasicLinalg::inner(m_dk, m_gkp1, size);
		return beta1 + beta2;
	}
}

} // namespace LNOT

#endif // LNOT_NON_LINEAR_CONJUGATE_GRADIENT_IMPL_HPP
