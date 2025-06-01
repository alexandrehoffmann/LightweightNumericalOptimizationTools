#ifndef LNOT_COUPLED_LANCZOS_SOLVER_IMPL_HPP
#define LNOT_COUPLED_LANCZOS_SOLVER_IMPL_HPP

namespace LNOT
{
	
template<typename T>
void CoupledLanczosSolver<T>::clearWorkSpace()
{
	if (m_v != nullptr) { delete[] m_v; m_v = nullptr; }
	if (m_q != nullptr) { delete[] m_q; m_q = nullptr; }
	if (m_w != nullptr) { delete[] m_w; m_w = nullptr; }
	Base::m_workCapacity = 0;
}

template<typename T> template<typename Op, bool solveInPlace> 
void CoupledLanczosSolver<T>::solve_impl(const Op& H, const Scalar* __restrict__ g, const Size size, std::bool_constant<solveInPlace>, Scalar* __restrict__ x) requires (IsHessianOp<Op>::value)
{
	if (Base::m_workCapacity < size)
	{
		clearWorkSpace();
		Base::m_workCapacity = size;
		m_v = new Scalar[Base::m_workCapacity];
		m_q = new Scalar[Base::m_workCapacity];
		m_w = new Scalar[Base::m_workCapacity];
	}
	if constexpr (solveInPlace)
	{
		H(x, m_v);
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_v[i] = -g[i] - m_v[i]; } // r = -b - Hx
	}
	else
	{
		std::fill(x, x+size, 0);
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_v[i] = -g[i]; } 
	}
	const Scalar normR0 = BasicLinalg::norm(m_v, size);
	
	Scalar beta_old = 0;
	Scalar l_old = 0;
	Scalar eta = -normR0;
	Scalar tau_old = 0;
	
	BasicLinalg::scal(1. / eta, size, m_v);
	
	std::copy(m_v, m_v + size, m_q);
	
	m_normR = normR0;
	
	const Scalar tol = Base::m_tol*m_normR;
	Base::m_info = Info::FAILURE;
	if (Base::m_out) { fmt::print(Base::m_out, "Coupled recurrence Lanczos solver : \n#Iteration residual tol\n"); }
	for (Base::m_nIt=0; Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e}\n", Base::m_nIt, m_normR, tol); }
		if (m_normR < tol) { Base::m_info = Info::SUCCESS; break; }
		H(m_q, m_w);
		
		const Scalar tau   = BasicLinalg::inner(m_q, m_w, size);
		const Scalar alpha = tau + tau_old*l_old*l_old;
		// solve T_k h_k = -\|r_0\|e_1
		const Scalar d = alpha - beta_old*l_old;
		if (d < std::numeric_limits<Scalar>::epsilon()) { Base::m_info = Info::NEGATIVE_CURVATURE; break; }
		const Scalar invD = 1. / d;
		BasicLinalg::axpy(invD*eta, m_q, size, x);
		// Resume Lanczos iteration
		BasicLinalg::axpy(-tau, m_v, size, m_w);
		const Scalar beta = BasicLinalg::norm(m_w, size);
		// prepare next solve
		const Scalar l = invD*beta;
		
		eta     *= -l;
		beta_old = beta;
		l_old    = l;
		tau_old  = tau;
		m_normR  = std::abs(eta); 
		const Scalar invBeta = 1. / beta;
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_v[i] = invBeta*m_w[i]; m_q[i] = m_v[i] - l*m_q[i]; } 
	}
}
	
} // namespace LNOT

#include <LNOT/LinearSolvers/CoupledLanczosSolver_impl.hpp>

#endif // LNOT_COUPLED_LANCZOS_SOLVER_IMPL_HPP
