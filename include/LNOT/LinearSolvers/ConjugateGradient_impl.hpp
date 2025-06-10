#ifndef LNOT_CONJUGATE_GRADIENT_IMPL_HPP
#define LNOT_CONJUGATE_GRADIENT_IMPL_HPP

#include <LNOT/LinearSolvers/ConjugateGradient.hpp>
#include <LNOT/BasicLinalg/SymmetricDenseMatrixOp.hpp>
#include <LNOT/BasicLinalg/DiagonalPreconditionerOp.hpp>
#include <LNOT/BasicLinalg/BasicLinalg.hpp>

#include <fmt/core.h>
#include <fmt/format.h>

namespace LNOT
{

//// explicit template instanciations ////

extern template class ConjugateGradient<float>;
extern template class ConjugateGradient<double>;

// explicit instanciation for float
extern template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::true_type, Scalar*);

extern template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::false_type, Scalar*);

extern template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::true_type, Scalar*);

extern template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void ConjugateGradient<float>::solve_impl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<float>&, const Scalar*, const Size, std::false_type, Scalar*);

// explicit instanciation for double
extern template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::true_type, Scalar*);

extern template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const Scalar*, const Size, std::false_type, Scalar*);

extern template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);
extern template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::true_type, Scalar*);

extern template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);
extern template void ConjugateGradient<double>::solve_impl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&, const DiagonalPreconditionerOp<double>&, const Scalar*, const Size, std::false_type, Scalar*);

//// method implementations ////

template<typename T> 
void ConjugateGradient<T>::clearWorkSpace()
{
	if (m_z  != nullptr) { delete[] m_z;  m_z  = nullptr; }
	if (m_r  != nullptr) { delete[] m_r;  m_r  = nullptr; }
	if (m_p  != nullptr) { delete[] m_p;  m_p  = nullptr; }
	if (m_Hp != nullptr) { delete[] m_Hp; m_Hp = nullptr; }
	Base::m_workCapacity = 0;
}

template<typename T> 
void ConjugateGradient<T>::resizeWorkSpace(const Size newSize)
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

template<typename T> template<typename Op, bool solveInPlace> 
void ConjugateGradient<T>::solve_impl(const Op& H, const Scalar* __restrict__ g, const Size size, std::bool_constant<solveInPlace>, Scalar* __restrict__ x) requires (IsHessianOp<Op>::value)
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
	std::copy(m_r, m_r + size, m_p);
	
	m_sqNormR = BasicLinalg::squaredNorm(m_r, size);
	
	const Scalar tol2 = Base::getSquaredResidualThreshold();
	
	Base::m_info = Info::FAILURE;
	if (Base::m_out) { fmt::print(Base::m_out, "#CG solver : \n#Iteration residual tol\n"); }
	for (Base::m_nIt=0; Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e}\n", Base::m_nIt, std::sqrt(m_sqNormR), std::sqrt(tol2)); }
		if (m_sqNormR < tol2) { Base::m_info = Info::SUCCESS; break; }
		H(m_p, m_Hp);
		const Scalar alpha = m_sqNormR / BasicLinalg::inner(m_p, m_Hp, size);
		
		if (alpha < std::numeric_limits<Scalar>::epsilon()) { Base::m_info = Info::NEGATIVE_CURVATURE; break; }
		
		BasicLinalg::axpy( alpha, m_p,  size,   x);
		BasicLinalg::axpy(-alpha, m_Hp, size, m_r);
		
		const Scalar sqNormNewR = BasicLinalg::squaredNorm(m_r, size);
		const Scalar beta = sqNormNewR / m_sqNormR;
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_p[i] = m_r[i] + beta*m_p[i]; } 
		
		m_sqNormR = sqNormNewR;
	}
}

template<typename T> template<typename HesOp, typename PrecOp, bool solveInPlace> 
void ConjugateGradient<T>::solve_impl(const HesOp& H, const PrecOp& invB, const Scalar* g, const Size size, std::bool_constant<solveInPlace>, Scalar* x) requires (AreHessianOps<HesOp,PrecOp>::value)
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
	
	m_sqNormR = BasicLinalg::inner(m_r, m_z, size);
	
	const Scalar tol2 = Base::getSquaredResidualThreshold();
	
	Base::m_info = Info::FAILURE;
	if (Base::m_out) { fmt::print(Base::m_out, "#Preconditioned CG solver : \n#Iteration residual tol\n"); }
	for (Base::m_nIt=0; Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e}\n", Base::m_nIt, std::sqrt(m_sqNormR), std::sqrt(tol2)); }
		if (m_sqNormR < tol2) { Base::m_info = Info::SUCCESS; break; }
		H(m_p, m_Hp);
		const Scalar alpha = m_sqNormR / BasicLinalg::inner(m_p, m_Hp, size);
		
		if (alpha < std::numeric_limits<Scalar>::epsilon()) { Base::m_info = Info::NEGATIVE_CURVATURE; break; }
		
		BasicLinalg::axpy( alpha, m_p,  size,   x);
		BasicLinalg::axpy(-alpha, m_Hp, size, m_r);
		
		invB(m_r, m_z);
		
		const Scalar sqNormNewR = BasicLinalg::inner(m_r, m_z, size);
		const Scalar beta = sqNormNewR / m_sqNormR;
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_p[i] = m_z[i] + beta*m_p[i]; } 
		
		m_sqNormR = sqNormNewR;
	}
}

} // namespace LNOT

#endif // LNOT_CONJUGATE_GRADIENT_IMPL_HPP
