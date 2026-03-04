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

// explicit instanciation for float
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);

extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);

extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);

extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);

// explicit instanciation for double
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);

extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);

extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);

extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);

// explicit instanciation for double
extern template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);
extern template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);
extern template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);
extern template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);

extern template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);

extern template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);
extern template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);
extern template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);
extern template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,true>, Scalar*);

extern template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);
extern template void LanczosSolver<long double>::solveImpl(const SymmetricDenseMatrixOp<long double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<long double>&, const Scalar*, const Size, BIC::Fixed<bool,false>, Scalar*);

//// method implementations ////

template<typename T> 
void LanczosSolver<T>::clearWorkSpace()
{
	if (m_Bv_old != nullptr) { delete[] m_Bv_old; m_Bv_old = nullptr; }
	if (m_Bv     != nullptr) { delete[] m_Bv;     m_Bv     = nullptr; }
	if (m_v      != nullptr) { delete[] m_v;      m_v      = nullptr; }
	if (m_p      != nullptr) { delete[] m_p;      m_p      = nullptr; }
	if (m_w      != nullptr) { delete[] m_w;      m_w      = nullptr; }
	m_workCapacity = 0;
}

template<typename T> 
void LanczosSolver<T>::resizeWorkSpace(const Size newSize)
{
	if (m_workCapacity < newSize)
	{
		clearWorkSpace();
		m_workCapacity = newSize;
		m_Bv_old  = new Scalar[m_workCapacity];
		m_Bv      = new Scalar[m_workCapacity];
		m_v       = new Scalar[m_workCapacity];
		m_p       = new Scalar[m_workCapacity];
		m_w  = new Scalar[m_workCapacity];
	}
}

template<typename T> template<typename HesOp, typename PrecOp, typename ASize, typename Bool>
void LanczosSolver<T>::solveImpl(const HesOp& H, const PrecOp& invB, const Scalar* g, const ASize size, Bool solveInPlace, Scalar* x) requires(isHessianOp<HesOp> and isHessianOp<PrecOp> and isSize<ASize>)
{
	resizeWorkSpace(size);
	if constexpr (solveInPlace)
	{
		H(x, m_Bv);
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_Bv[i] = -g[i] - m_Bv[i]; } // r = -b - Hx
	}
	else
	{
		std::fill(x, x+size, 0);
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_Bv[i] = -g[i]; } 
	}
	std::fill(m_Bv_old, m_Bv_old + size, 0);
	invB(m_Bv, m_v);
	
	const Scalar precNormR0 = sqrt(BasicLinalg::inner(m_Bv, m_v, size));
	
	const FPComparator<Scalar> cmp;
	
	Scalar beta_old{};
	Scalar l_old{};
	Scalar eta = -precNormR0;
	
	BasicLinalg::scal(Scalar(1) / eta, size, m_Bv);
	BasicLinalg::scal(Scalar(1) / eta, size, m_v);
	
	std::fill(m_p, m_p + size, 0);
	
	m_precNormR = precNormR0;
	
	const Scalar relTol = m_relTol*m_precNormR;
	
	m_info = Info::FAILURE;
	if (m_out) { fmt::print(m_out, "#Preconditioned Lanczos solver : \n#Iteration residual relative_tol absolute_tol\n"); }
	for (m_nIt=0; m_nIt!=m_maxIt; ++m_nIt)
	{
		if (m_out) { fmt::print(m_out, "{} {:10.2e} {:10.2e} {:10.2e}\n", m_nIt, m_precNormR, relTol, m_absTol); }
		if (m_precNormR < relTol or m_precNormR < m_absTol) { m_info = Info::SUCCESS; break; }
		
		H(m_v, m_w); 
		// m_w = Hv_{k}
		const Scalar alpha = BasicLinalg::inner(m_v, m_w, size);
		// solve T_k h_k = -\|r_0\|e_1
		const Scalar d = alpha - beta_old*l_old;
		if (not cmp.isDefPositive(d)) { m_info = Info::NEGATIVE_CURVATURE; break; }
		const Scalar invD = Scalar(1) / d;
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_p[i] = invD*(m_v[i] - beta_old*m_p[i]); } 
		BasicLinalg::axpy(eta, m_p, size, x);
		// resume Lanczos iteration
		BasicLinalg::axpbypz(-alpha, m_Bv, -beta_old, m_Bv_old, size, m_w);
		// m_w = \beta_{k}Bv_{k+1}
		invB(m_w, m_v);
		const Scalar beta    = sqrt(BasicLinalg::inner(m_w, m_v, size));
		const Scalar invBeta = Scalar(1) / beta;
		BasicLinalg::scal(invBeta, size, m_v);
		std::copy(m_Bv, m_Bv + size, m_Bv_old);
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
