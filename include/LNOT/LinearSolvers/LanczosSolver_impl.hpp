#ifndef LNOT_LANCZOS_SOLVER_IMPL_HPP
#define LNOT_LANCZOS_SOLVER_IMPL_HPP

#include <LNOT/LinearSolvers/LanczosSolver.hpp>
#include <LNOT/BasicLinalg/SymmetricDenseMatrixOp.hpp>
#include <LNOT/BasicLinalg/DiagonalPreconditionerOp.hpp>
#include <LNOT/BasicLinalg/BasicLinalg.hpp>
#include <LNOT/FloatingPoint/FPComparator.hpp>

#include <fmt/core.h>
#include <fmt/format.h>

namespace LNOT
{

//// explicit template instanciations ////

extern template class LanczosSolver<float>;
extern template class LanczosSolver<double>;
extern template class LanczosSolver<long double>;

//// method implementations ////

template<typename T> 
void LanczosSolver<T>::resizeWorkSpace(const Size newSize)
{
	if (not m_Bv_old or m_workCapacity < newSize)
	{
		m_workCapacity = newSize;
		m_Bv_old = std::make_unique<Scalar[]>(this->m_workCapacity);
		m_Bv     = std::make_unique<Scalar[]>(this->m_workCapacity);
		m_v      = std::make_unique<Scalar[]>(this->m_workCapacity);
		m_p      = std::make_unique<Scalar[]>(this->m_workCapacity);
		m_w      = std::make_unique<Scalar[]>(this->m_workCapacity);
	}
}

template<typename T> template<typename HesOp, typename PrecOp, typename ASize, typename Bool>
void LanczosSolver<T>::solveImpl(const HesOp& H, const PrecOp& invB, const Scalar* g, const ASize size, Bool solveInPlace, Scalar* x) requires(isHessianOp<HesOp> and isHessianOp<PrecOp> and isSize<ASize>)
{
	resizeWorkSpace(size);
	if constexpr (solveInPlace)
	{
		H(x, m_Bv.get());
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_Bv[i] = -g[i] - m_Bv[i]; } // r = -b - Hx
	}
	else
	{
		std::fill(x, x+size, 0);
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_Bv[i] = -g[i]; } 
	}
	std::fill(m_Bv_old.get(), m_Bv_old.get() + size, 0);
	invB(m_Bv.get(), m_v.get());
	
	const Scalar precNormR0 = sqrt(BasicLinalg::inner(m_Bv.get(), m_v.get(), size));
	
	constexpr FPComparator<Scalar> cmp;
	
	Scalar beta_old{};
	Scalar l_old{};
	Scalar eta = -precNormR0;
	
	BasicLinalg::scal(Scalar(1) / eta, size, m_Bv.get());
	BasicLinalg::scal(Scalar(1) / eta, size, m_v.get());
	
	std::fill(m_p.get(), m_p.get() + size, 0);
	
	m_precNormR = precNormR0;
	
	const Scalar relTol = m_relTol*m_precNormR;
	
	m_info = Info::FAILURE;
	if (m_out) { fmt::println(m_out, "#Preconditioned Lanczos solver : \n#Iteration residual relative_tol absolute_tol"); }
	for (m_nIt=0; m_nIt!=m_maxIt; ++m_nIt)
	{
		if (m_out) { fmt::println(m_out, "{} {:10.2e} {:10.2e} {:10.2e}", m_nIt, m_precNormR, relTol, m_absTol); std::fflush(m_out); }
		if (m_precNormR < relTol or m_precNormR < m_absTol) { m_info = Info::SUCCESS; break; }
		
		H(m_v.get(), m_w.get()); 
		// m_w = Hv_{k}
		const Scalar alpha = BasicLinalg::inner(m_v.get(), m_w.get(), size);
		// solve T_k h_k = -\|r_0\|e_1
		const Scalar d = alpha - beta_old*l_old;
		if (not cmp.isDefPositive(d)) { m_info = Info::NEGATIVE_CURVATURE; break; }
		const Scalar invD = Scalar(1) / d;
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_p[i] = invD*(m_v[i] - beta_old*m_p[i]); } 
		BasicLinalg::axpy(eta, m_p.get(), size, x);
		// resume Lanczos iteration
		BasicLinalg::axpbypz(-alpha, m_Bv.get(), -beta_old, m_Bv_old.get(), size, m_w.get());
		// m_w = \beta_{k}Bv_{k+1}
		invB(m_w.get(), m_v.get());
		const Scalar beta    = sqrt(BasicLinalg::inner(m_w.get(), m_v.get(), size));
		const Scalar invBeta = Scalar(1) / beta;
		BasicLinalg::scal(invBeta, size, m_v.get());
		std::copy(m_Bv.get(), m_Bv.get() + size, m_Bv_old.get());
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_Bv[i] = invBeta*m_w[i]; }
		// prepare next solve
		const Scalar l = invD*beta;
		
		eta     *= -l;
		beta_old = beta;
		l_old    = l;
		m_precNormR  = std::abs(eta); 
	}
}
	
} // namespace LNOT


#endif // LNOT_LANCZOS_SOLVER_IMPL_HPP
