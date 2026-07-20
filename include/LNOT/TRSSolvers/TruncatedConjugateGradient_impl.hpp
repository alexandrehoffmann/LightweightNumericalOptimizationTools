#ifndef LNOT_TRUNCATED_CONJUGATE_GRADIENT_IMPL_HPP
#define LNOT_TRUNCATED_CONJUGATE_GRADIENT_IMPL_HPP

#include <LNOT/TRSSolvers/TruncatedConjugateGradient.hpp>
#include <LNOT/BasicLinalg/BasicLinalg.hpp>
#include <LNOT/FloatingPoint/FPComparator.hpp>

#include <fmt/core.h>
#include <fmt/format.h>

namespace LNOT
{
	
//// explicit template instanciations ////

extern template class TruncatedConjugateGradient<float>;
extern template class TruncatedConjugateGradient<double>;

//// method implementations ////
	
template<typename T> 
void TruncatedConjugateGradient<T>::resizeWorkSpace(const Size newSize)
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

template<typename T> template<typename HesOp, typename PrecOp, typename ASize> 
auto TruncatedConjugateGradient<T>::solveImpl(HesOp&& H, PrecOp&& invB, const Scalar* g, const ASize size, const Scalar& delta, Scalar* x) -> Scalar requires(areHessianOps<HesOp,PrecOp> and isSize<ASize>)
{	
	using AdlMath::sqrt;
	
	resizeWorkSpace(size);
	
	std::fill(x, x+size, 0);
	m_modelReduction = 0;
	Scalar sqNormX{};
	
	#pragma omp simd
	for (Size i=0; i!=size; ++i) { m_r[i] = -g[i]; } 
	invB.eval(m_r.get(), m_z.get());
	
	std::copy(m_z.get(), m_z.get() + size, m_p.get());
	
	m_precSqNormR = BasicLinalg::inner(m_r.get(), m_z.get(), size);
	
	Scalar precSqNormP = m_precSqNormR;
	Scalar precInnerXP{}; // x = 0, thus (x, p)_B = 0
	
	const Scalar relTol2 = m_relTol*m_relTol*m_precSqNormR;
	const Scalar absTol2 = m_absTol*m_absTol;
	const Scalar delta2  = delta*delta;
	
	constexpr FPComparator<Scalar> cmp;
	const     FPComparator<Scalar> cmpTr(m_relTolTr, m_absTolTr);
	
	m_info = Info::FAILURE;
	if (m_out) { fmt::println(m_out, "#Truncated Preconditioned CG solver : \n#Iteration residual relative_tol absolute_tol"); }
	for (m_nIt=0; m_nIt!=m_maxIt; ++m_nIt)
	{
		if (m_out) { fmt::println(m_out, "{} {:10.2e} {:10.2e} {:10.2e}", m_nIt, sqrt(m_precSqNormR), sqrt(relTol2), sqrt(m_absTol)); std::fflush(m_out); }
 		if (m_precSqNormR < relTol2 or m_precSqNormR < absTol2) { m_info = Info::SUCCESS; break; }
 		
 		H.eval(m_p.get(), m_Hp.get());
 		
 		const Scalar alpha = m_precSqNormR / BasicLinalg::inner(m_p.get(), m_Hp.get(), size);
 		const Scalar sqNormXnext = sqNormX + 2*alpha*precInnerXP + alpha*alpha*precSqNormP;
 		
 		if ((not cmp.isDefPositive(alpha)) or cmpTr.isDefGreaterThan(sqNormXnext, delta2)) 
 		{
			// we need to find tau > 0 such that |x + tau*p|^2 = delta^2
			const Scalar tau = getPolyMaxRoot(precSqNormP, 2*precInnerXP, sqNormX - delta*delta);
			
			sqNormX = delta;
			
			if (cmp.isDefNegative(tau)) { m_info = Info::BREAKDOWN; break; } 
			
			m_modelReduction = BasicLinalg::updateModelReduction(m_modelReduction, x, g, tau, m_p.get(), m_Hp.get(), size);
			BasicLinalg::axpy(tau, m_p.get(), size, x);
			
			break;
		}
		
 		m_modelReduction = BasicLinalg::updateModelReduction(m_modelReduction, x, g, alpha, m_p.get(), m_Hp.get(), size);
		
		BasicLinalg::axpy( alpha, m_p.get(),  size, x);
		BasicLinalg::axpy(-alpha, m_Hp.get(), size, m_r.get());
		
		invB.eval(m_r.get(), m_z.get());
		
		const Scalar precSqNormNewR = BasicLinalg::inner(m_r.get(), m_z.get(), size);
		const Scalar beta = precSqNormNewR / m_precSqNormR;
		
		sqNormX     = sqNormXnext;
		precInnerXP = beta*precInnerXP + BasicLinalg::inner(x, m_r.get(), size) + alpha*beta*precSqNormP;
		precSqNormP = precSqNormNewR + beta*BasicLinalg::inner(m_p.get(), m_r.get(), size) + beta*beta*precSqNormP;
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_p[i] = m_z[i] + beta*m_p[i]; } 
		
		m_precSqNormR = precSqNormNewR;
	}
	
	return sqrt(sqNormX);
}

template<typename T>
auto TruncatedConjugateGradient<T>::getPolyMaxRoot(const Scalar a, const Scalar b, const Scalar c) -> Scalar
{ // solve for x > 0 a*x*x + b*x + c
	using AdlMath::sqrt;
	
	const Scalar delta = b*b - 4*a*c;
	if (delta < NumTraits<Scalar>::epsilon) { return -b / (2*a); }
	const Scalar x1 = (-b + sqrt(delta)) / (2*a);
	const Scalar x2 = (-b - sqrt(delta)) / (2*a);
	return std::max(x1, x2);
}
	
} // namespace LNOT

#endif // LNOT_TRUNCATED_CONJUGATE_GRADIENT_IMPL_HPP
