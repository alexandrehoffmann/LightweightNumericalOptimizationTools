#ifndef LNOT_BFGS_IMPL_HPP
#define LNOT_BFGS_IMPL_HPP

#include <LNOT/GeneralSolvers/BFGS.hpp>
#include <LNOT/Containers/CircularBuffer.hpp>
#include <LNOT/LineSearches/BisectionLineSearch.hpp>
#include <LNOT/LineSearches/BacktrackingLineSearch.hpp>

#include <span>
#include <fmt/ranges.h>

namespace LNOT
{

//// explicit template instanciations ////

extern template class BFGS< BisectionLineSearch<float> >;
extern template class BFGS< BacktrackingLineSearch<float> >;

extern template class BFGS< BisectionLineSearch<double> >;
extern template class BFGS< BacktrackingLineSearch<double> >;

extern template class BFGS< BisectionLineSearch<long double> >;
extern template class BFGS< BacktrackingLineSearch<long double> >;

//// method implementations ////

template<typename LineSearch> 
void BFGS<LineSearch>::clearWorkSpaceImpl()
{
	if (m_gk    != nullptr) { delete[] m_gk;    m_gk    = nullptr; }
	if (m_gkp1  != nullptr) { delete[] m_gkp1;  m_gkp1  = nullptr; }
	if (m_sk    != nullptr) { delete[] m_sk;    m_sk    = nullptr; }
	if (m_yk    != nullptr) { delete[] m_yk;    m_yk    = nullptr; }
	if (m_uk    != nullptr) { delete[] m_uk;    m_uk    = nullptr; }
	if (m_invBk != nullptr) { delete[] m_invBk; m_invBk = nullptr; }
	Base::m_workCapacity = 0;
}

template<typename LineSearch> template<FirstOrderOracle_concept Oracle, bool solveInPlace> 
void BFGS<LineSearch>::solveImpl(Oracle& oracle, std::bool_constant<solveInPlace>, Scalar* x)
{
	using Oracle_Size = typename Oracle::Size;
	
	const Oracle_Size size = oracle.getNDims();
	
	if (Base::m_workCapacity < size)
	{
		clearWorkSpaceImpl();
		Base::m_workCapacity = size;

		m_gk    = new Scalar[Base::m_workCapacity];
		m_gkp1  = new Scalar[Base::m_workCapacity];
		m_sk    = new Scalar[Base::m_workCapacity];
		m_yk    = new Scalar[Base::m_workCapacity];
		m_uk    = new Scalar[Base::m_workCapacity];
		m_invBk = new Scalar[Base::m_workCapacity*Base::m_workCapacity];
	}
	if constexpr (not solveInPlace) { std::fill(x, x + size, 0); }
	
	std::fill(m_invBk, m_invBk + size*size, 0);
	for (Size i=0; i!=size; ++i) { m_invBk[i + i*size] = 1; }
	
	Base::m_innerIts.clear();
	
	oracle.setCurrentPoint(x);
	oracle.getGradient(m_gk);
	
	Base::m_fx = oracle.getValue();
	Base::m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk, size);
	
	const Scalar tol2 = Base::m_tol*Base::m_tol*std::max(Scalar(1), Base::m_squaredNormGrad);
	
	if (Base::m_out != nullptr) { fmt::print(Base::m_out, "#BFGS method\n#Iteration f(x) residual tol\n"); }
	
	Base::m_info = Info::FAILURE;
	for (Base::m_nIt=0;Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e} {:10.2e}\n", Base::m_nIt, Base::m_fx, Base::m_squaredNormGrad, tol2); }
		if (Base::m_squaredNormGrad < tol2) { Base::m_info = Info::SUCCESS; break; }
		
		BasicLinalg::symMatrixVectorProd(StorageOrder::ROW_MAJOR, UpLo::LOWER, Scalar(-1), m_invBk, m_gk, size, std::false_type{}, m_sk);
		Base::m_innerIts.push_back(1);
		
		////
		const Scalar step_len = m_lineSearch.solve(x, Base::m_fx, m_gk, m_sk, oracle);
		
		if (step_len < std::numeric_limits<Scalar>::epsilon()) { Base::m_info = Info::BREAKDOWN; break; }
		
		BasicLinalg::scal(step_len, size, m_sk);
		BasicLinalg::axpy(Scalar(1), m_sk, size, x);
		
		oracle.setCurrentPoint(x);
		oracle.getGradient(m_gkp1);
		
		#pragma omp simd
		for (Size j=0; j!=size; ++j) 
		{ 
			m_yk[j] = m_gkp1[j] - m_gk[j];
		}
		
		const Scalar sk_dot_yk = BasicLinalg::inner(m_sk, m_yk, size);
		
		if (sk_dot_yk > 0)
		{			
			BasicLinalg::symMatrixVectorProd(StorageOrder::ROW_MAJOR, UpLo::LOWER, Scalar(1), m_invBk, m_yk, size, std::false_type{}, m_uk);
			
			const Scalar yk_dot_uk = BasicLinalg::inner(m_yk, m_uk, size);
			
			BasicLinalg::symRk1Update(StorageOrder::ROW_MAJOR, UpLo::LOWER, (sk_dot_yk + yk_dot_uk) / (sk_dot_yk*sk_dot_yk), m_sk, size, m_invBk);
			BasicLinalg::symRk2Update(StorageOrder::ROW_MAJOR, UpLo::LOWER, Scalar(-1) / sk_dot_yk, m_uk, m_sk, size, m_invBk);
		}
		
		////
		
		std::copy(m_gkp1, m_gkp1 + size, m_gk);
		Base::m_fx = oracle.getValue();
		Base::m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk, size);
	}
}

} // namespace LNOT

#endif // LNOT_BFGS_IMPL_HPP
