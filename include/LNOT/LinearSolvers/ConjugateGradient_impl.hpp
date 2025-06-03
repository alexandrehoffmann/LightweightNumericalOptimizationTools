#ifndef LNOT_CONJUGATE_GRADIENT_IMPL_HPP
#define LNOT_CONJUGATE_GRADIENT_IMPL_HPP

#include <LNOT/LinearSolvers/ConjugateGradient.hpp>
#include <LNOT/BasicLinalg.hpp>

#include <fmt/core.h>
#include <fmt/format.h>

namespace LNOT
{

template<typename T> 
void ConjugateGradient<T>::clearWorkSpace()
{
	if (m_r  != nullptr) { delete[] m_r;  m_r  = nullptr; }
	if (m_p  != nullptr) { delete[] m_p;  m_p  = nullptr; }
	if (m_Hp != nullptr) { delete[] m_Hp; m_Hp = nullptr; }
	Base::m_workCapacity = 0;
}

template<typename T> template<typename Op, bool solveInPlace> 
void ConjugateGradient<T>::solve_impl(const Op& H, const Scalar* __restrict__ g, const Size size, std::bool_constant<solveInPlace>, Scalar* __restrict__ x) requires (IsHessianOp<Op>::value)
{
	if (Base::m_workCapacity < size)
	{
		clearWorkSpace();
		Base::m_workCapacity = size;
		m_r  = new Scalar[Base::m_workCapacity];
		m_p  = new Scalar[Base::m_workCapacity];
		m_Hp = new Scalar[Base::m_workCapacity];
	}
	
	if constexpr (solveInPlace)
	{
		H(x, m_r);
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_r[i] = -g[i] - m_r[i]; } // r = -b - Hx
	}
	else
	{
		std::fill(x, x+size, 0);
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_r[i] = -g[i]; } 
	}
	std::copy(m_r, m_r + size, m_p);
	
	m_sqNormR = BasicLinalg::squaredNorm(m_r, size);
	
	const Scalar tol2 = Base::m_tol*Base::m_tol*m_sqNormR;
	
	Base::m_info = Info::FAILURE;
	if (Base::m_out) { fmt::print(Base::m_out, "#CG solver : \n#Iteration residual tol\n"); }
	for (Base::m_nIt=0; Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e}\n", Base::m_nIt, std::sqrt(m_sqNormR), std::sqrt(tol2)); }
		if (m_sqNormR < tol2) { Base::m_info = Info::SUCCESS; break; }
		H(m_p, m_Hp);
		const Scalar alpha = m_sqNormR / BasicLinalg::inner(m_p, m_Hp, size);
		
		if (alpha < std::numeric_limits<Scalar>::epsilon()) { Base::m_info = Info::NEGATIVE_CURVATURE; break; }
		
		BasicLinalg::axpy( alpha, m_p,  size,   x);
		BasicLinalg::axpy(-alpha, m_Hp, size, m_r);
		
		const Scalar sqNormNewR = BasicLinalg::squaredNorm(m_r, size);
		const Scalar beta = sqNormNewR / m_sqNormR;
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_p[i] = m_r[i] + beta*m_p[i]; } 
		
		m_sqNormR = sqNormNewR;
	}
}

} // namespace LNOT

#endif // LNOT_CONJUGATE_GRADIENT_IMPL_HPP
