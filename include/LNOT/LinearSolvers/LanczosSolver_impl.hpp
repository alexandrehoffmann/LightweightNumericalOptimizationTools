#ifndef LNOT_LANCZOS_SOLVER_IMPL_HPP
#define LNOT_LANCZOS_SOLVER_IMPL_HPP

#include <LNOT/LinearSolvers/LanczosSolver.hpp>
#include <LNOT/BasicLinalg/SymmetricDenseMatrixOp.hpp>
#include <LNOT/BasicLinalg/DiagonalPreconditionerOp.hpp>
#include <LNOT/BasicLinalg/BasicLinalg.hpp>

#include <fmt/core.h>
#include <fmt/format.h>

namespace LNOT
{

//// explicit template instanciations ////

extern template class LanczosSolver<float>;
extern template class LanczosSolver<double>;

// explicit instanciation for float
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);

extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);

extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);

extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void LanczosSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);

// explicit instanciation for double
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);

extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const IdentityPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);

extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);

extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void LanczosSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);

//// method implementations ////

template<typename T> 
void LanczosSolver<T>::clearWorkSpace()
{
	if (m_Bv_old != nullptr) { delete[] m_Bv_old; m_Bv_old = nullptr; }
	if (m_Bv     != nullptr) { delete[] m_Bv;     m_Bv     = nullptr; }
	if (m_v      != nullptr) { delete[] m_v;      m_v      = nullptr; }
	if (m_p      != nullptr) { delete[] m_p;      m_p      = nullptr; }
	if (m_w      != nullptr) { delete[] m_w;      m_w      = nullptr; }
	Base::m_workCapacity = 0;
}

template<typename T> 
void LanczosSolver<T>::resizeWorkSpace(const Size newSize)
{
	if (Base::m_workCapacity < newSize)
	{
		clearWorkSpace();
		Base::m_workCapacity = newSize;
		m_Bv_old  = new Scalar[Base::m_workCapacity];
		m_Bv      = new Scalar[Base::m_workCapacity];
		m_v       = new Scalar[Base::m_workCapacity];
		m_p       = new Scalar[Base::m_workCapacity];
		m_w  = new Scalar[Base::m_workCapacity];
	}
}

template<typename T> template<typename HesOp, typename PrecOp, bool solveInPlace> 
void LanczosSolver<T>::solveImpl(const HesOp& H, const PrecOp& invB, const Scalar* g, const Size size, std::bool_constant<solveInPlace>, Scalar* x) requires (AreHessianOps<HesOp,PrecOp>::value)
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
	
	const Scalar precNormR0 = std::sqrt(BasicLinalg::inner(m_Bv, m_v, size));
	
	Scalar beta_old = 0;
	Scalar l_old = 0;
	Scalar eta = -precNormR0;
	
	BasicLinalg::scal(Scalar(1) / eta, size, m_Bv);
	BasicLinalg::scal(Scalar(1) / eta, size, m_v);
	
	std::fill(m_p, m_p + size, 0);
	
	m_precNormR = precNormR0;
	
	const Scalar tol = Base::getResidualThreshold();
	
	Base::m_info = Info::FAILURE;
	if (Base::m_out) { fmt::print(Base::m_out, "#Preconditioned Lanczos solver : \n#Iteration residual tol\n"); }
	for (Base::m_nIt=0; Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e}\n", Base::m_nIt, m_precNormR, tol); }
		if (m_precNormR < tol) { Base::m_info = Info::SUCCESS; break; }
		
		H(m_v, m_w); 
		// m_w = Hv_{k}
		const Scalar alpha = BasicLinalg::inner(m_v, m_w, size);
		// solve T_k h_k = -\|r_0\|e_1
		const Scalar d = alpha - beta_old*l_old;
		if (d < std::numeric_limits<Scalar>::epsilon()) { Base::m_info = Info::NEGATIVE_CURVATURE; break; }
		const Scalar invD = Scalar(1) / d;
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_p[i] = invD*(m_v[i] - beta_old*m_p[i]); } 
		BasicLinalg::axpy(eta, m_p, size, x);
		fmt::print("eta = {}\n", eta);
		// resume Lanczos iteration
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_w[i] += -alpha*m_Bv[i] - beta_old*m_Bv_old[i]; } 
		// m_w = \beta_{k}Bv_{k+1}
		invB(m_w, m_v);
		const Scalar beta    = std::sqrt(BasicLinalg::inner(m_w, m_v, size));
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
