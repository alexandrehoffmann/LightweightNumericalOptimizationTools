#ifndef LNOT_COUPLED_LANCZOS_SOLVER_IMPL_HPP
#define LNOT_COUPLED_LANCZOS_SOLVER_IMPL_HPP

#include <LNOT/LinearSolvers/CoupledLanczosSolver.hpp>
#include <LNOT/SymmetricDenseMatrixOp.hpp>
#include <LNOT/DiagonalPreconditionerOp.hpp>
#include <LNOT/BasicLinalg.hpp>

#include <fmt/core.h>
#include <fmt/format.h>

namespace LNOT
{

//// explicit template instanciations ////

extern template class CoupledLanczosSolver<float>;
extern template class CoupledLanczosSolver<double>;

// explicit instanciation for float
extern template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::true_type, Scalar*);

extern template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::false_type, Scalar*);

extern template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);

extern template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void CoupledLanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);

// explicit instanciation for double
extern template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::true_type, Scalar*);

extern template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::false_type, Scalar*);

extern template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);

extern template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void CoupledLanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);

//// method implementations ////

template<typename T>
void CoupledLanczosSolver<T>::clearWorkSpace()
{
	if (m_v     != nullptr) { delete[] m_v;     m_v     = nullptr; }
	if (m_Bv     != nullptr) { delete[] m_Bv;     m_Bv     = nullptr; }
	if (m_q     != nullptr) { delete[] m_q;     m_q     = nullptr; }
	if (m_hat_Bv != nullptr) { delete[] m_hat_Bv; m_hat_Bv = nullptr; }
	Base::m_workCapacity = 0;
}

template<typename T>
void CoupledLanczosSolver<T>::resizeWorkSpace(const Size newSize)
{
	if (Base::m_workCapacity < newSize)
	{
		clearWorkSpace();
		Base::m_workCapacity = newSize;
		m_v     = new Scalar[Base::m_workCapacity];
		m_Bv     = new Scalar[Base::m_workCapacity];
		m_q     = new Scalar[Base::m_workCapacity];
		m_hat_Bv = new Scalar[Base::m_workCapacity];
	}
}

template<typename T> template<typename Op, bool solveInPlace> 
void CoupledLanczosSolver<T>::solve_impl(const Op& H, const Scalar* __restrict__ g, const Size size, std::bool_constant<solveInPlace>, Scalar* __restrict__ x) requires (IsHessianOp<Op>::value)
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
	const Scalar normR0 = BasicLinalg::norm(m_Bv, size);
	
	Scalar beta_old = 0;
	Scalar l_old = 0;
	Scalar eta = -normR0;
	Scalar tau_old = 0;
	
	BasicLinalg::scal(Scalar(1) / eta, size, m_Bv);
	
	std::copy(m_Bv, m_Bv + size, m_q);
	
	m_normR = normR0;
	
	const Scalar tol = Base::getResidualThreshold();
	
	Base::m_info = Info::FAILURE;
	if (Base::m_out) { fmt::print(Base::m_out, "#Coupled recurrence Lanczos solver : \n#Iteration residual tol\n"); }
	for (Base::m_nIt=0; Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e}\n", Base::m_nIt, m_normR, tol); }
		if (m_normR < tol) { Base::m_info = Info::SUCCESS; break; }
		
		H(m_q, m_hat_Bv);
		const Scalar tau   = BasicLinalg::inner(m_q, m_hat_Bv, size);
		const Scalar alpha = tau + tau_old*l_old*l_old;
		// solve T_k h_k = -\|r_0\|e_1
		const Scalar d = alpha - beta_old*l_old;
		if (d < std::numeric_limits<Scalar>::epsilon()) { Base::m_info = Info::NEGATIVE_CURVATURE; break; }
		const Scalar invD = Scalar(1) / d;
		BasicLinalg::axpy(invD*eta, m_q, size, x);
		// Resume Lanczos iteration
		BasicLinalg::axpy(-tau, m_Bv, size, m_hat_Bv);
		const Scalar beta = BasicLinalg::norm(m_hat_Bv, size);
		// prepare next solve
		const Scalar l = invD*beta;
		
		eta     *= -l;
		beta_old = beta;
		l_old    = l;
		tau_old  = tau;
		m_normR  = std::abs(eta); 
		const Scalar invBeta = Scalar(1) / beta;
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_Bv[i] = invBeta*m_hat_Bv[i]; m_q[i] = m_Bv[i] - l*m_q[i]; } 
	}
}

template<typename T> template<typename HesOp, typename PrecOp, bool solveInPlace> 
void CoupledLanczosSolver<T>::solve_impl(const HesOp& H, const PrecOp& invB, const Scalar* __restrict__ g, const Size size, std::bool_constant<solveInPlace>, Scalar* __restrict__ x) requires (AreHessianOps<HesOp,PrecOp>::value)
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
	invB(m_Bv, m_v);
	
	const Scalar normR0 = std::sqrt(BasicLinalg::inner(m_Bv, m_v, size));	
	
	Scalar tau_old  = 0;
	Scalar beta_old = 0;
	Scalar l_old    = 0;
	Scalar eta      = -normR0;
	
	BasicLinalg::scal(Scalar(1) / eta, size, m_Bv);
	BasicLinalg::scal(Scalar(1) / eta, size, m_v);
	
	std::copy(m_v, m_v + size, m_q);
	m_normR = normR0;
	
	const Scalar tol = Base::getResidualThreshold();
	
	Base::m_info = Info::FAILURE;
	if (Base::m_out) { fmt::print(Base::m_out, "#Preconditioned Coupled recurrence Lanczos solver : \n#Iteration residual tol\n"); }
	for (Base::m_nIt=0; Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e}\n", Base::m_nIt, m_normR, tol); }
		if (m_normR < tol) { Base::m_info = Info::SUCCESS; break; }
		H(m_q, m_hat_Bv);
		const Scalar tau   = BasicLinalg::inner(m_q, m_hat_Bv, size);
		const Scalar alpha = tau + tau_old*l_old*l_old;
		// solve T_k h_k = -\|r_0\|e_1
		const Scalar d = alpha - beta_old*l_old;
		if (d < std::numeric_limits<Scalar>::epsilon()) { Base::m_info = Info::NEGATIVE_CURVATURE; break; }
		const Scalar invD = Scalar(1) / d;
		BasicLinalg::axpy(invD*eta, m_q, size, x);
		// Resume Lanczos iteration
		BasicLinalg::axpy(-tau, m_Bv, size, m_hat_Bv);
		invB(m_hat_Bv, m_v);
		const Scalar beta = std::sqrt(BasicLinalg::inner(m_hat_Bv, m_v, size));
		// prepare next solve
		const Scalar l = invD*beta;
		
		eta     *= -l;
		beta_old = beta;
		l_old    = l;
		tau_old  = tau;
		m_normR  = std::abs(eta); 
		const Scalar invBeta = Scalar(1) / beta;
		BasicLinalg::scal(invBeta, size, m_v);
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_Bv[i] = invBeta*m_hat_Bv[i]; m_q[i] = m_v[i] - l*m_q[i]; } 
	}
}
	
} // namespace LNOT

#include <LNOT/LinearSolvers/CoupledLanczosSolver_impl.hpp>

#endif // LNOT_COUPLED_LANCZOS_SOLVER_IMPL_HPP
