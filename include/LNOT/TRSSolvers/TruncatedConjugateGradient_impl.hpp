#ifndef LNOT_TRUNCATED_CONJUGATE_GRADIENT_IMPL_HPP
#define LNOT_TRUNCATED_CONJUGATE_GRADIENT_IMPL_HPP

#include <LNOT/TRSSolvers/TruncatedConjugateGradient.hpp>
#include <LNOT/BasicLinalg/SymmetricDenseMatrixOp.hpp>
#include <LNOT/BasicLinalg/BasicLinalg.hpp>

#include <fmt/core.h>
#include <fmt/format.h>

namespace LNOT
{
	
//// explicit template instanciations ////

extern template class TruncatedConjugateGradient<float>;
extern template class TruncatedConjugateGradient<double>;

// explicit instanciation for float
extern template void TruncatedConjugateGradient<float>::solve(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void TruncatedConjugateGradient<float>::solve(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void TruncatedConjugateGradient<float>::solve(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void TruncatedConjugateGradient<float>::solve(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);

// explicit instanciation for double
extern template void TruncatedConjugateGradient<double>::solve(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void TruncatedConjugateGradient<double>::solve(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void TruncatedConjugateGradient<double>::solve(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void TruncatedConjugateGradient<double>::solve(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
	
//// method implementations ////
	
template<typename T>
void TruncatedConjugateGradient<T>::clearWorkSpace()
{
	if (m_r  != nullptr) { delete[] m_r;  } 
	if (m_p  != nullptr) { delete[] m_p;  } 
	if (m_Hp != nullptr) { delete[] m_Hp; }
	Base::m_workCapacity = 0;
}

template<typename T> template<typename Op> 
void TruncatedConjugateGradient<T>::solve(const Op& H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) requires (IsHessianOp<Op>::value)
{
	if (Base::m_workCapacity < size)
	{
		clearWorkSpace();
		Base::m_workCapacity = size;
		m_r  = new Scalar[Base::m_workCapacity];
		m_p  = new Scalar[Base::m_workCapacity];
		m_Hp = new Scalar[Base::m_workCapacity];
	}
	
	std::fill(x, x+size, 0);
	Base::m_modelReduction = 0;
	Scalar sqNormX = 0;
	#pragma omp simd
	for (Size i=0; i!=size; ++i) { m_r[i] = -g[i]; } 
	
	std::copy(m_r, m_r + size, m_p);
	
	m_sqNormR = BasicLinalg::squaredNorm(m_r, size);
	
	const Scalar tol2 = Base::m_tol*Base::m_tol*m_sqNormR;
	const Scalar deltaTol2 = (delta + Base::m_tolTr)*(delta + Base::m_tolTr);
	
	Base::m_info = Info::FAILURE;
	if (Base::m_out) { fmt::print(Base::m_out, "#Truncated CG solver : \n#Iteration residual tol\n"); }
	for (Base::m_nIt=0; Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e}\n", Base::m_nIt, std::sqrt(m_sqNormR), std::sqrt(tol2)); }
		if (m_sqNormR < tol2) { Base::m_info = Info::SUCCESS; break; }
		
		H(m_p, m_Hp);
		
		const Scalar alpha = m_sqNormR / BasicLinalg::inner(m_p, m_Hp, size);
		const Scalar sqNormXnext = sqNormX + 2*alpha*BasicLinalg::inner(x, m_p, size) + alpha*alpha*BasicLinalg::squaredNorm(m_p, size);
		
		if (alpha < std::numeric_limits<Scalar>::epsilon() or sqNormXnext >= deltaTol2) 
		{ 
			// we need to find tau > 0 such that |x + tau*p|^2 = delta^2
			const Scalar tau = getPolyMaxRoot(BasicLinalg::squaredNorm(m_p, size), 2*BasicLinalg::inner(x, m_p, size), sqNormX - delta*delta);

			if (tau < 0) { Base::m_info = Info::BREAKDOWN; break; } 
			
			#pragma omp simd reduction(+:Base::m_modelReduction)
			for (Size i=0;i !=size; ++i) { Base::m_modelReduction += tau*(x[i]*m_Hp[i] + Scalar(0.5)*tau*m_p[i]*m_Hp[i] + g[i]*m_p[i]); }
			BasicLinalg::axpy(tau, m_p, size, x);
			break;
		}
		else
		{
			#pragma omp simd reduction(+:Base::m_modelReduction)
			for (Size i=0;i !=size; ++i) { Base::m_modelReduction += alpha*(x[i]*m_Hp[i] + Scalar(0.5)*alpha*m_p[i]*m_Hp[i] + g[i]*m_p[i]); }
		}
		
		BasicLinalg::axpy( alpha, m_p,  size, x);
		BasicLinalg::axpy(-alpha, m_Hp, size, m_r);
		
		sqNormX = sqNormXnext;
		
		const Scalar sqNormNewR = BasicLinalg::squaredNorm(m_r, size);
		const Scalar beta = sqNormNewR / m_sqNormR;
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_p[i] = m_r[i] + beta*m_p[i]; } 
		
		m_sqNormR = sqNormNewR;
	}
}

template<typename T>
auto TruncatedConjugateGradient<T>::getPolyMaxRoot(const Scalar a, const Scalar b, const Scalar c) -> Scalar
{ // solve for x > 0 a*x*x + b*x + c
	const Scalar delta = b*b - 4*a*c;
	if (delta < std::numeric_limits<Scalar>::epsilon()) { return -b / (2*a); }
	const Scalar x1 = (-b + std::sqrt(delta)) / (2*a);
	const Scalar x2 = (-b - std::sqrt(delta)) / (2*a);
	return std::max(x1, x2);
}
	
} // namespace LNOT

#endif // LNOT_TRUNCATED_CONJUGATE_GRADIENT_IMPL_HPP
