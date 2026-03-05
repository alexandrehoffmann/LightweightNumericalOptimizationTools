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

// explicit instanciation for float
extern template void ConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);
extern template void ConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);
extern template void ConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);
extern template void ConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);

extern template void ConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);
extern template void ConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);
extern template void ConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);
extern template void ConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);

extern template void ConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);
extern template void ConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);
extern template void ConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);
extern template void ConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);

extern template void ConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);
extern template void ConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);
extern template void ConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);
extern template void ConjugateGradient<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);

// explicit instanciation for double
extern template void ConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);
extern template void ConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);
extern template void ConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);
extern template void ConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);

extern template void ConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);
extern template void ConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);
extern template void ConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);
extern template void ConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);

extern template void ConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);
extern template void ConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);
extern template void ConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);
extern template void ConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);

extern template void ConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);
extern template void ConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);
extern template void ConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);
extern template void ConjugateGradient<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);

// explicit instanciation for long double
extern template void ConjugateGradient<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);
extern template void ConjugateGradient<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);
extern template void ConjugateGradient<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);
extern template void ConjugateGradient<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);

extern template void ConjugateGradient<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);
extern template void ConjugateGradient<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);
extern template void ConjugateGradient<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);
extern template void ConjugateGradient<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, const BIC::Fixed<bool,false>, Scalar*);

extern template void ConjugateGradient<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);
extern template void ConjugateGradient<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);
extern template void ConjugateGradient<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);
extern template void ConjugateGradient<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, const BIC::Fixed<bool,true>, Scalar*);

extern template void ConjugateGradient<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void ConjugateGradient<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void ConjugateGradient<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void ConjugateGradient<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);

//// method implementations ////

template<typename T> 
void ConjugateGradient<T>::clearWorkSpace()
{
	if (m_z  != nullptr) { delete[] m_z;  m_z  = nullptr; }
	if (m_r  != nullptr) { delete[] m_r;  m_r  = nullptr; }
	if (m_p  != nullptr) { delete[] m_p;  m_p  = nullptr; }
	if (m_Hp != nullptr) { delete[] m_Hp; m_Hp = nullptr; }
	m_workCapacity = 0;
}

template<typename T> 
void ConjugateGradient<T>::resizeWorkSpace(const Size newSize)
{
	if (m_workCapacity < newSize)
	{
		clearWorkSpace();
		m_workCapacity = newSize;
		m_z  = new Scalar[m_workCapacity];
		m_r  = new Scalar[m_workCapacity];
		m_p  = new Scalar[m_workCapacity];
		m_Hp = new Scalar[m_workCapacity];
	}
}

template<typename T> template<typename HesOp, typename PrecOp, typename ASize, typename Bool>
void ConjugateGradient<T>::solveImpl(const HesOp& H, const PrecOp& invB, const Scalar* g, const ASize size, const Bool solveInPlace, Scalar* x) requires (isHessianOp<HesOp> and isHessianOp<PrecOp> and isSize<ASize>)
{		
	resizeWorkSpace(size);
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
	invB(m_r, m_z);
	std::copy(m_z, m_z + size, m_p);
	
	m_precSqNormR = BasicLinalg::inner(m_r, m_z, size);
	
	const Scalar relTol2 = m_relTol*m_relTol*m_precSqNormR;
	const Scalar absTol2 = m_absTol*m_absTol;
	
	const FPComparator<Scalar> cmp;
	
	m_info = Info::FAILURE;
	if (m_out) { fmt::println(m_out, "#Preconditioned CG solver : \n#Iteration residual relative_tol absolute_tol"); }
	for (m_nIt=0; m_nIt!=m_maxIt; ++m_nIt)
	{
		if (m_out) { fmt::println(m_out, "{} {:10.2e} {:10.2e} {:10.2e}", m_nIt, sqrt(m_precSqNormR), sqrt(relTol2), sqrt(absTol2)); std::fflush(m_out); }
		if (m_precSqNormR < relTol2 or m_precSqNormR < absTol2) { m_info = Info::SUCCESS; break; }
		H(m_p, m_Hp);
		const Scalar alpha = m_precSqNormR / BasicLinalg::inner(m_p, m_Hp, size);
		
		if (not cmp.isDefPositive(alpha)) { m_info = Info::NEGATIVE_CURVATURE; break; }
		
		BasicLinalg::axpy( alpha, m_p,  size,   x);
		BasicLinalg::axpy(-alpha, m_Hp, size, m_r);
		
		invB(m_r, m_z);
		
		const Scalar precSqNormNewR = BasicLinalg::inner(m_r, m_z, size);
		const Scalar beta = precSqNormNewR / m_precSqNormR;
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_p[i] = m_z[i] + beta*m_p[i]; } 
		
		m_precSqNormR = precSqNormNewR;
	}
}

} // namespace LNOT

#endif // LNOT_CONJUGATE_GRADIENT_IMPL_HPP
