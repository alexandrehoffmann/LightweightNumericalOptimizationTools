#ifndef LNOT_LANCZOS_SOLVER_IMPL_HPP
#define LNOT_LANCZOS_SOLVER_IMPL_HPP

#include <LNOT/LinearSolvers/LanczosSolver.hpp>
#include <LNOT/SymmetricDenseMatrixOp.hpp>
#include <LNOT/BasicLinalg.hpp>

#include <fmt/core.h>
#include <fmt/format.h>

namespace LNOT
{

//// explicit template instanciations ////

extern template class LanczosSolver<float>;
extern template class LanczosSolver<double>;

// explicit instanciation for float
extern template void LanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);
extern template void LanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);
extern template void LanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);
extern template void LanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);

extern template void LanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);
extern template void LanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);
extern template void LanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);
extern template void LanczosSolver<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);

// explicit instanciation for double
extern template void LanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);
extern template void LanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);
extern template void LanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);
extern template void LanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::true_type, Scalar* x);

extern template void LanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);
extern template void LanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);
extern template void LanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);
extern template void LanczosSolver<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const Scalar* g, const Size size, std::false_type, Scalar* x);

//// method implementations ////

template<typename T> 
void LanczosSolver<T>::clearWorkSpace()
{
	if (m_v_old != nullptr) { delete[] m_v_old; m_v_old = nullptr; }
	if (m_v     != nullptr) { delete[] m_v;     m_v     = nullptr; }
	if (m_p     != nullptr) { delete[] m_p;     m_p     = nullptr; }
	if (m_w     != nullptr) { delete[] m_w;     m_w     = nullptr; }
	Base::m_workCapacity = 0;
}

template<typename T>  template<typename Op, bool solveInPlace> 
void LanczosSolver<T>::solve_impl(const Op& H, const Scalar* __restrict__ g, const Size size, std::bool_constant<solveInPlace>, Scalar* __restrict__ x) requires (IsHessianOp<Op>::value)
{
	if (Base::m_workCapacity < size)
	{
		clearWorkSpace();
		Base::m_workCapacity = size;
		m_v_old  = new Scalar[Base::m_workCapacity];
		m_v      = new Scalar[Base::m_workCapacity];
		m_p      = new Scalar[Base::m_workCapacity];
		m_w      = new Scalar[Base::m_workCapacity];
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
	
	BasicLinalg::scal(Scalar(1) / eta, size, m_v);
	
	std::fill(m_v_old, m_v_old + size, 0);
	std::fill(m_p,     m_p     + size, 0);
	
	m_normR = normR0;
	
	const Scalar tol = Base::m_tol*m_normR;
	
	Base::m_info = Info::FAILURE;
	if (Base::m_out) { fmt::print(Base::m_out, "#Lanczos solver : \n#Iteration residual tol\n"); }
	for (Base::m_nIt=0; Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e}\n", Base::m_nIt, m_normR, tol); }
		if (m_normR < tol) { Base::m_info = Info::SUCCESS; break; }
		
		H(m_v, m_w);
		const Scalar alpha = BasicLinalg::inner(m_v, m_w, size);
		// solve T_k h_k = -\|r_0\|e_1
		const Scalar d = alpha - beta_old*l_old;
		if (d < std::numeric_limits<Scalar>::epsilon()) { Base::m_info = Info::NEGATIVE_CURVATURE; break; }
		const Scalar invD = Scalar(1) / d;
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_p[i] = invD*(m_v[i] - beta_old*m_p[i]); } 
		BasicLinalg::axpy(eta, m_p, size, x);
		// resume Lanczos iteration
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_w[i] += -alpha*m_v[i] - beta_old*m_v_old[i]; } 
		const Scalar beta = BasicLinalg::norm(m_w, size);
		// prepare next solve
		const Scalar l = invD*beta;
		
		eta     *= -l;
		beta_old = beta;
		l_old    = l;
		m_normR  = std::abs(eta); 
		const Scalar invBeta = Scalar(1) / beta;
		#pragma omp simd
		for (Size i=0; i!=size; ++i) 
		{ 
			m_v_old[i] = m_v[i];
			m_v[i]     = invBeta*m_w[i];
		} 
	}
}
	
} // namespace LNOT


#endif // LNOT_LANCZOS_SOLVER_IMPL_HPP
