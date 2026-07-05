#ifndef LNOT_CONJUGATE_GRADIENT_IMPL_HPP
#define LNOT_CONJUGATE_GRADIENT_IMPL_HPP

#include <LNOT/LinearSolvers/ConjugateGradient.hpp>
#include <LNOT/BasicLinalg/SymmetricDenseMatrixOp.hpp>
#include <LNOT/BasicLinalg/DiagonalPreconditionerOp.hpp>
#include <LNOT/BasicLinalg/BasicLinalg.hpp>
#include <LNOT/FloatingPoint/FPComparator.hpp>

#include <fmt/core.h>
#include <fmt/format.h>

namespace LNOT
{

//// explicit template instanciations ////

extern template class ConjugateGradient<float>;
extern template class ConjugateGradient<double>;
extern template class ConjugateGradient<long double>;

//// method implementations ////

template<typename T> 
void ConjugateGradient<T>::resizeWorkSpace(const Size newSize)
{
	if (not m_z or m_workCapacity < newSize)
	{
		m_workCapacity = newSize;
		m_z  = std::make_unique<Scalar[]>(m_workCapacity);
		m_r  = std::make_unique<Scalar[]>(m_workCapacity);
		m_p  = std::make_unique<Scalar[]>(m_workCapacity);
		m_Hp = std::make_unique<Scalar[]>(m_workCapacity);
	}
}

template<typename T> template<typename HesOp, typename PrecOp, typename ASize, typename Bool>
void ConjugateGradient<T>::solveImpl(const HesOp& H, const PrecOp& invB, const Scalar* g, const ASize size, const Bool solveInPlace, Scalar* x) requires (isHessianOp<HesOp> and isHessianOp<PrecOp> and isSize<ASize>)
{		
	resizeWorkSpace(size);
	if constexpr (solveInPlace)
	{
		H(x, m_r.get());
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_r[i] = -g[i] - m_r[i]; } // r = -b - Hx
	}
	else
	{
		std::fill(x, x+size, 0);
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_r[i] = -g[i]; } 
	}
	invB(m_r.get(), m_z.get());
	std::copy(m_z.get(), m_z.get() + size, m_p.get());
	
	m_precSqNormR = BasicLinalg::inner(m_r.get(), m_z.get(), size);
	
	const Scalar relTol2 = m_relTol*m_relTol*m_precSqNormR;
	const Scalar absTol2 = m_absTol*m_absTol;
	
	constexpr FPComparator<Scalar> cmp;
	
	m_info = Info::FAILURE;
	if (m_out) { fmt::println(m_out, "#Preconditioned CG solver : \n#Iteration residual relative_tol absolute_tol"); }
	for (m_nIt=0; m_nIt!=m_maxIt; ++m_nIt)
	{
		if (m_out) { fmt::println(m_out, "{} {:10.2e} {:10.2e} {:10.2e}", m_nIt, sqrt(m_precSqNormR), sqrt(relTol2), sqrt(absTol2)); std::fflush(m_out); }
		if (m_precSqNormR < relTol2 or m_precSqNormR < absTol2) { m_info = Info::SUCCESS; break; }
		H(m_p.get(), m_Hp.get());
		const Scalar alpha = m_precSqNormR / BasicLinalg::inner(m_p.get(), m_Hp.get(), size);
		
		if (not cmp.isDefPositive(alpha)) { m_info = Info::NEGATIVE_CURVATURE; break; }
		
		BasicLinalg::axpy( alpha, m_p.get(),  size,   x);
		BasicLinalg::axpy(-alpha, m_Hp.get(), size, m_r.get());
		
		invB(m_r.get(), m_z.get());
		
		Scalar precSqNormNewR = BasicLinalg::inner(m_r.get(), m_z.get(), size);
		
		const Scalar beta = precSqNormNewR / m_precSqNormR;
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_p[i] = m_z[i] + beta*m_p[i]; } 
		
		m_precSqNormR = std::move(precSqNormNewR);
	}
}

} // namespace LNOT

#endif // LNOT_CONJUGATE_GRADIENT_IMPL_HPP
