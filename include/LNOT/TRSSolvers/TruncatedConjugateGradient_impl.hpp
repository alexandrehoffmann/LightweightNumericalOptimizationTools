#ifndef LNOT_TRUNCATED_CONJUGATE_GRADIENT_IMPL_HPP
#define LNOT_TRUNCATED_CONJUGATE_GRADIENT_IMPL_HPP

#include <LNOT/TRSSolvers/TruncatedConjugateGradient.hpp>
#include <LNOT/BasicLinalg/SymmetricDenseMatrixOp.hpp>
#include <LNOT/BasicLinalg/BasicLinalg.hpp>
#include <LNOT/BasicLinalg/DiagonalPreconditionerOp.hpp>

#include <fmt/core.h>
#include <fmt/format.h>

namespace LNOT
{
	
//// explicit template instanciations ////

extern template class TruncatedConjugateGradient<float>;
extern template class TruncatedConjugateGradient<double>;

// explicit instanciation for float
extern template void TruncatedConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void TruncatedConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void TruncatedConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void TruncatedConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);

extern template void TruncatedConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void TruncatedConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void TruncatedConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void TruncatedConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);

// explicit instanciation for double
extern template void TruncatedConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void TruncatedConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void TruncatedConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void TruncatedConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);

extern template void TruncatedConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void TruncatedConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void TruncatedConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void TruncatedConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);

//// method implementations ////
	
template<typename T>
void TruncatedConjugateGradient<T>::clearWorkSpace()
{
	if (m_z  != nullptr) { delete[] m_z;  } 
	if (m_r  != nullptr) { delete[] m_r;  } 
	if (m_p  != nullptr) { delete[] m_p;  } 
	if (m_Hp != nullptr) { delete[] m_Hp; }
	Base::m_workCapacity = 0;
}

template<typename T> 
void TruncatedConjugateGradient<T>::resizeWorkSpace(const Size newSize)
{
	if (Base::m_workCapacity < newSize)
	{
		clearWorkSpace();
		Base::m_workCapacity = newSize;
		m_z  = new Scalar[Base::m_workCapacity];
		m_r  = new Scalar[Base::m_workCapacity];
		m_p  = new Scalar[Base::m_workCapacity];
		m_Hp = new Scalar[Base::m_workCapacity];
	}
}

template<typename T> template<typename HesOp, typename PrecOp, typename ASize> 
void TruncatedConjugateGradient<T>::solveImpl(const HesOp& H, const PrecOp& invB, const Scalar* g, const ASize size, const Scalar& delta, Scalar* x) requires(AreHessianOps<HesOp,PrecOp>::value and IsSize<ASize>::value)
{	
	resizeWorkSpace(size);
	
	std::fill(x, x+size, 0);
	Base::m_modelReduction = 0;
	Scalar sqNormX = 0;
	
	#pragma omp simd
	for (Size i=0; i!=size; ++i) { m_r[i] = -g[i]; } 
	invB(m_r, m_z);
	
	std::copy(m_z, m_z + size, m_p);
	
	m_precSqNormR = BasicLinalg::inner(m_r, m_z, size);
	
	Scalar precSqNormP = m_precSqNormR;
	Scalar precInnerXP = 0; // x = 0, thus (x, p)_B = 0
	
	const Scalar tol2      = Base::getSquaredResidualThreshold();
	const Scalar deltaTol2 = (delta + Base::m_tolTr)*(delta + Base::m_tolTr);
	
	Base::m_info = Info::FAILURE;
	if (Base::m_out) { fmt::print(Base::m_out, "#Truncated Preconditioned CG solver : \n#Iteration residual tol\n"); }
	for (Base::m_nIt=0; Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e}\n", Base::m_nIt, std::sqrt(m_precSqNormR), std::sqrt(tol2)); }
 		if (m_precSqNormR < tol2) { Base::m_info = Info::SUCCESS; break; }
 		
 		H(m_p, m_Hp);
 		
 		const Scalar alpha = m_precSqNormR / BasicLinalg::inner(m_p, m_Hp, size);
 		const Scalar sqNormXnext = sqNormX + 2*alpha*precInnerXP + alpha*alpha*precSqNormP;
 		
 		if (alpha < std::numeric_limits<Scalar>::epsilon() or sqNormXnext >= deltaTol2) 
 		{
			// we need to find tau > 0 such that |x + tau*p|^2 = delta^2
			const Scalar tau = getPolyMaxRoot(precSqNormP, 2*precInnerXP, sqNormX - delta*delta);
			
			if (tau < 0) { Base::m_info = Info::BREAKDOWN; break; } 
			
			for (Size i=0;i !=size; ++i) { Base::m_modelReduction += tau*(x[i]*m_Hp[i] + Scalar(0.5)*tau*m_p[i]*m_Hp[i] + g[i]*m_p[i]); }
			BasicLinalg::axpy(tau, m_p, size, x);
			break;
		}
		
 		for (Size i=0;i !=size; ++i) { Base::m_modelReduction += alpha*(x[i]*m_Hp[i] + Scalar(0.5)*alpha*m_p[i]*m_Hp[i] + g[i]*m_p[i]); }
		
		BasicLinalg::axpy( alpha, m_p,  size, x);
		BasicLinalg::axpy(-alpha, m_Hp, size, m_r);
		
		invB(m_r, m_z);
		
		const Scalar precSqNormNewR = BasicLinalg::inner(m_r, m_z, size);
		const Scalar beta = precSqNormNewR / m_precSqNormR;
		
		sqNormX     = sqNormXnext;
		precInnerXP = beta*precInnerXP + BasicLinalg::inner(x, m_r, size) + alpha*beta*precSqNormP;
		precSqNormP = precSqNormNewR + beta*BasicLinalg::inner(m_p, m_r, size) + beta*beta*precSqNormP;
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_p[i] = m_z[i] + beta*m_p[i]; } 
		
		m_precSqNormR = precSqNormNewR;
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
