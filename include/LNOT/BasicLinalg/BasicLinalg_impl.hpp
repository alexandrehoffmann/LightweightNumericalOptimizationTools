#ifndef LNOT_BASIC_LINALG_IMPL_HPP
#define LNOT_BASIC_LINALG_IMPL_HPP

#include <LNOT/BasicLinalg/BasicLinalg.hpp>
#include <BIC/Core.hpp>

#include <numeric>

namespace LNOT::BasicLinalg
{

template<typename Scalar, typename Size> 
void stridedCopy(const Scalar* x, const Size xStride, const Size N, Scalar* y, const Size yStride)
{
#ifdef LNOT_WITH_BLAS
	if constexpr      (std::is_same<Scalar, float>::value)  { return cblas_scopy(blasint(N), x, blasint(xStride), y, blasint(yStride)); }
	else if constexpr (std::is_same<Scalar, double>::value) { return cblas_dcopy(blasint(N), x, blasint(xStride), y, blasint(yStride)); }
	else
	{
#endif // LNOT_WITH_BLAS
		#pragma omp simd 
		for (BIC::Mutable<Size> i=0; i!=N; ++i) 
		{ 
			y[i*yStride] = x[i*xStride]; 
		} 
#ifdef LNOT_WITH_BLAS
	}
#endif // LNOT_WITH_BLAS	
}

template<typename Scalar, typename Size>
Scalar squaredNorm(const Scalar* x, const Size N) 
{ 
#ifdef LNOT_WITH_BLAS
	if constexpr      (std::is_same<Scalar, float>::value)  { return cblas_sdot(blasint(N), x, 1, x, 1); }
	else if constexpr (std::is_same<Scalar, double>::value) { return cblas_ddot(blasint(N), x, 1, x, 1); }
	else
	{
#endif // LNOT_WITH_BLAS
		Scalar res(0); 
		for (BIC::Mutable<Size> i=0; i!=N; ++i) { res += x[i]*x[i]; } 
		return res; 
#ifdef LNOT_WITH_BLAS
	}
#endif // LNOT_WITH_BLAS
}

template<typename Scalar, typename Size> 
Scalar norm(const Scalar* x, const Size N) 
{
#ifdef LNOT_WITH_BLAS
	if constexpr      (std::is_same<Scalar, float>::value)  { return cblas_snrm2(blasint(N), x, 1); }
	else if constexpr (std::is_same<Scalar, double>::value) { return cblas_dnrm2(blasint(N), x, 1); }
	else
	{
#endif // LNOT_WITH_BLAS
		return sqrt(squaredNorm(x, N));
#ifdef LNOT_WITH_BLAS
	}
#endif // LNOT_WITH_BLAS
} 

template<typename Scalar, typename Size>
Scalar weightedSquaredNorm(const Scalar* x, const Scalar* w, const Size N)
{
	Scalar res(0);
	
	for (BIC::Mutable<Size> i=0; i!=N; ++i) {  res += x[i]*x[i]*w[i]; }
	
	return res;
}

template<typename Scalar, typename Size>
Scalar inner(const Scalar* x, const Scalar* y, const Size N)
{ 
#ifdef LNOT_WITH_BLAS
	if constexpr      (std::is_same<Scalar, float>::value)  { return cblas_sdot(blasint(N), x, 1, y, 1); }
	else if constexpr (std::is_same<Scalar, double>::value) { return cblas_ddot(blasint(N), x, 1, y, 1); }
	else
	{
#endif // LNOT_WITH_BLAS
		Scalar res(0); 
		for (BIC::Mutable<Size> i=0; i!=N; ++i) { res += x[i]*y[i]; } 
		return res; 
#ifdef LNOT_WITH_BLAS
	}
#endif // LNOT_WITH_BLAS
}

template<typename Scalar, typename Size>
Scalar weightedInner(const Scalar* x, const Scalar* y, const Scalar* w, const Size N)
{ 
	Scalar res(0); 
	for (BIC::Mutable<Size> i=0; i!=N; ++i) { res += x[i]*y[i]*w[i]; } 
	return res; 
}

template<typename Scalar, typename Size>
void axpy(const Scalar alpha, const Scalar* x, const Size N, Scalar* y)
{
#ifdef LNOT_WITH_BLAS
	if constexpr      (std::is_same<Scalar, float>::value)  { cblas_saxpy(blasint(N), alpha, x, 1, y, 1); }
	else if constexpr (std::is_same<Scalar, double>::value) { cblas_daxpy(blasint(N), alpha, x, 1, y, 1); }
	else
	{
#else
#endif // LNOT_WITH_BLAS	
		#pragma omp simd 
		for (BIC::Mutable<Size> i=0; i!=N; ++i)
		{
			y[i] += alpha*x[i];
		}
#ifdef LNOT_WITH_BLAS
	}
#endif // LNOT_WITH_BLAS
}

template<typename Scalar, typename Size>
void scal(const Scalar alpha, const Size N, Scalar* x)
{
#ifdef LNOT_WITH_BLAS
	if constexpr      (std::is_same<Scalar, float>::value)  { cblas_sscal(blasint(N), alpha, x, 1); }
	else if constexpr (std::is_same<Scalar, double>::value) { cblas_dscal(blasint(N), alpha, x, 1); }
	else
	{
#endif // LNOT_WITH_BLAS
		#pragma omp simd 
		for (BIC::Mutable<Size> i=0; i!=N; ++i)
		{
			x[i] *= alpha;
		}
#ifdef LNOT_WITH_BLAS
	}
#endif // LNOT_WITH_BLAS
}

template<typename Scalar, typename Size, typename Bool> 
void symMatrixVectorProd(const StorageOrder layout, const UpLo uplo, const Scalar alpha, const Scalar* A, const Scalar* x, const Size N, const Bool incrY, Scalar* y)
{
#ifdef LNOT_WITH_BLAS
	if      constexpr (std::is_same<Scalar, float>::value)  { cblas_ssymv(CBLAS_ORDER(layout), CBLAS_UPLO(uplo), blasint(N), alpha, A, blasint(N), x, 1, blasint(incrY), y, 1); }
	else if constexpr (std::is_same<Scalar, double>::value) { cblas_dsymv(CBLAS_ORDER(layout), CBLAS_UPLO(uplo), blasint(N), alpha, A, blasint(N), x, 1, blasint(incrY), y, 1); }
	else
	{
#endif // LNOT_WITH_BLAS
		if (not incrY) { std::fill(y, y + N, 0); }
		
		if      constexpr (std::is_same<Scalar, float>::value)       { lnot_symMatrixVectorProd_f (lnot_mat_StorageOrder(layout), lnot_mat_UpLo(uplo), alpha, A, x, lnot_Size(N), y); }
		else if constexpr (std::is_same<Scalar, double>::value)      { lnot_symMatrixVectorProd_d (lnot_mat_StorageOrder(layout), lnot_mat_UpLo(uplo), alpha, A, x, lnot_Size(N), y); }
		else if constexpr (std::is_same<Scalar, long double>::value) { lnot_symMatrixVectorProd_ld(lnot_mat_StorageOrder(layout), lnot_mat_UpLo(uplo), alpha, A, x, lnot_Size(N), y); }
		else
		{
			const Size iStride = layout == StorageOrder::ROW_MAJOR and uplo == UpLo::LOWER ? N : 1;
			const Size jStride = layout == StorageOrder::ROW_MAJOR and uplo == UpLo::LOWER ? 1 : N;
			
			for (BIC::Mutable<Size> i=0; i!=N; ++i)
			{
				for (BIC::Mutable<Size> j=0; j!=(i+1); ++j)
				{
					y[i] += alpha*A[i*iStride + j*jStride]*x[j];
				}
				for (BIC::Mutable<Size> j=i+1; j!=N; ++j)
				{
					y[i] += alpha*A[j*iStride + i*jStride]*x[j];
				}
			}
		}
#ifdef LNOT_WITH_BLAS
	}
#endif // LNOT_WITH_BLAS
}

template<typename Scalar, typename Size>
void symRk1Update(const StorageOrder layout, const UpLo uplo, const Scalar alpha, const Scalar* x, const Size N, Scalar* A)
{
#ifdef LNOT_WITH_BLAS
	if      constexpr (std::is_same<Scalar, float>::value)  { cblas_ssyr(CBLAS_ORDER(layout), CBLAS_UPLO(uplo), blasint(N), alpha, x, 1, A, blasint(N)); }
	else if constexpr (std::is_same<Scalar, double>::value) { cblas_dsyr(CBLAS_ORDER(layout), CBLAS_UPLO(uplo), blasint(N), alpha, x, 1, A, blasint(N)); }
#else
	if      constexpr (std::is_same<Scalar, float>::value)  { lnot_symRk1Update_f(lnot_mat_StorageOrder(layout), lnot_mat_UpLo(uplo), alpha, x, lnot_Size(N), A); }
	else if constexpr (std::is_same<Scalar, double>::value) { lnot_symRk1Update_d(lnot_mat_StorageOrder(layout), lnot_mat_UpLo(uplo), alpha, x, lnot_Size(N), A); }
#endif // LNOT_WITH_BLAS
	else if constexpr (std::is_same<Scalar, long double>::value) { lnot_symRk1Update_ld(lnot_mat_StorageOrder(layout), lnot_mat_UpLo(uplo), alpha, x, lnot_Size(N), A); }
	else
	{
		const Size iStride = layout == StorageOrder::ROW_MAJOR and uplo == UpLo::LOWER ? N : 1;
		const Size jStride = layout == StorageOrder::ROW_MAJOR and uplo == UpLo::LOWER ? 1 : N; 
	
		for (BIC::Mutable<Size> i=0; i!=N; ++i)
		{
			for (BIC::Mutable<Size> j=0; j!=(i+1); ++j)
			{
				A[i*iStride + j*jStride] += alpha*x[i]*x[j];
			}
		}
	}
}

template<typename Scalar, typename Size>
void symRk2Update(StorageOrder layout, UpLo uplo, const Scalar alpha, const Scalar* x, const Scalar* y, const Size N, Scalar* A)
{
#ifdef LNOT_WITH_BLAS
	if constexpr      (std::is_same<Scalar, float>::value)  { cblas_ssyr2(CBLAS_ORDER(layout), CBLAS_UPLO(uplo), blasint(N), alpha, x, 1, y, 1, A, blasint(N)); }
	else if constexpr (std::is_same<Scalar, double>::value) { cblas_dsyr2(CBLAS_ORDER(layout), CBLAS_UPLO(uplo), blasint(N), alpha, x, 1, y, 1, A, blasint(N)); }
#else
	if constexpr      (std::is_same<Scalar, float>::value)  { lnot_symRk2Update_f(lnot_mat_StorageOrder(layout), lnot_mat_UpLo(uplo), alpha, x, y, lnot_Size(N), A); }
	else if constexpr (std::is_same<Scalar, double>::value) { lnot_symRk2Update_d(lnot_mat_StorageOrder(layout), lnot_mat_UpLo(uplo), alpha, x, y, lnot_Size(N), A); }
#endif // LNOT_WITH_BLAS
	else if constexpr (std::is_same<Scalar, long double>::value) { lnot_symRk2Update_ld(lnot_mat_StorageOrder(layout), lnot_mat_UpLo(uplo), alpha, x, y, lnot_Size(N), A); }
	else
	{
		const Size iStride = layout == StorageOrder::ROW_MAJOR and uplo == UpLo::LOWER ? N : 1;
		const Size jStride = layout == StorageOrder::ROW_MAJOR and uplo == UpLo::LOWER ? 1 : N; 
	
		for (BIC::Mutable<Size> i=0; i!=N; ++i)
		{
			for (BIC::Mutable<Size> j=0; j!=(i+1); ++j)
			{
				A[i*iStride + j*jStride] += alpha*x[i]*y[j] + alpha*y[i]*x[j];
			}
		}
	}
}

namespace Tridiag
{

template<typename Scalar, typename Size>
Scalar norm1(const Scalar* alpha, const Scalar* beta, const Size N)
{
	if constexpr      (std::is_same<Scalar, float>::value)       { return lnot_tridiag_norm1_f (alpha, beta, lnot_Size(N)); }
	else if constexpr (std::is_same<Scalar, double>::value)      { return lnot_tridiag_norm1_d (alpha, beta, lnot_Size(N)); }
	else if constexpr (std::is_same<Scalar, long double>::value) { return lnot_tridiag_norm1_ld(alpha, beta, lnot_Size(N)); }
	else
	{
		if (N == 1) { return std::abs(alpha[0]); }

		Scalar res = std::abs(alpha[0]) + std::abs(beta[0]);
		#pragma omp simd reduction(max:res)
		for (BIC::Mutable<Size> i=1; i!=BIC::Mutable<Size>(N-1); ++i)
		{
			res = std::max(res, std::abs(alpha[i]) + std::abs(beta[i]) + std::abs(beta[i-1]));
		}
		return std::max(res, std::abs(alpha[N-1]) + std::abs(beta[N-2]));
	}
}

namespace LDLt
{

template<typename Scalar, typename Size>
bool compute(const Scalar* alpha, const Scalar* beta, const Size size, const Scalar shift, Scalar* invDelta, Scalar* l)
{
	constexpr Scalar epsilon = std::numeric_limits<Scalar>::epsilon();
	
	if constexpr (std::is_same<Scalar, float>::value)       { return lnot_tridiag_ldlt_compute_f (alpha, beta, size, shift, epsilon, invDelta, l); }
	if constexpr (std::is_same<Scalar, double>::value)      { return lnot_tridiag_ldlt_compute_d (alpha, beta, size, shift, epsilon, invDelta, l); }
	if constexpr (std::is_same<Scalar, long double>::value) { return lnot_tridiag_ldlt_compute_ld(alpha, beta, size, shift, epsilon, invDelta, l); }
	else
	{
		if ((alpha[0] + shift) < epsilon) { return false; }
		
		invDelta[0] = Scalar(1) / (alpha[0] + shift);
		if (size == 1) { return true; } 
		l[0] = beta[0]*invDelta[0];
		
		for (BIC::Mutable<Size> i=1; i!=BIC::Mutable<Size>(size-1);++i)
		{
			const Scalar delta_i = alpha[i] + shift - beta[i-1]*l[i-1];
			if (delta_i < epsilon) { return false; }
			invDelta[i] = Scalar(1) / delta_i;
			l[i] = beta[i]*invDelta[i];
		}
		{ 
			const Size i = Size(size-1);
			const Scalar delta_i = alpha[i] + shift - beta[i-1]*l[i-1];
			if (delta_i < epsilon) { return false; }
			invDelta[i] = Scalar(1) / delta_i;
		}
		
		return true;
	}
}

template<typename Scalar, typename Size>
void solveLowerUnit(const Scalar* l, const Size size, const Scalar b1, Scalar* x)
{
	x[0] = b1;
	for (BIC::Mutable<Size> i=1; i!=size; ++i) { x[i] = -l[i-1]*x[i-1]; }
}

template<typename Scalar, typename Size>
void solveInplaceLower(const Scalar* l, const Size size, Scalar* x)
{
	for (BIC::Mutable<Size> i=1; i!=size; ++i) { x[i] -= l[i-1]*x[i-1]; }
}

template<typename Scalar, typename Size>
void solveInplaceUpper(const Scalar* l, const Size size, Scalar* x)
{
	for (BIC::Mutable<Size> i=BIC::Mutable<Size>(size-2); i!=BIC::Mutable<Size>(-1); --i) { x[i] -= l[i]*x[i+1]; }
}

template<typename Scalar, typename Size>
void solveUnit(const Scalar* invD, const Scalar* l, const Size size, const Scalar b1, Scalar* x)
{
	if constexpr      (std::is_same<Scalar, float>::value)       { lnot_tridiag_ldlt_solveUnit_f (invD, l, lnot_Size(size), b1, x); }
	else if constexpr (std::is_same<Scalar, double>::value)      { lnot_tridiag_ldlt_solveUnit_d (invD, l, lnot_Size(size), b1, x); }
	else if constexpr (std::is_same<Scalar, long double>::value) { lnot_tridiag_ldlt_solveUnit_ld(invD, l, lnot_Size(size), b1, x); }
	else
	{
		solveLowerUnit(l, size, b1, x);     // Solve Lz = b
		#pragma omp simd
		for (BIC::Mutable<Size> i=0; i!=size; ++i) { x[i] *= invD[i]; }
		solveInplaceUpper(l, size, x);   // Solve L^T x = y
	}
}

} // namespace LDLt

} // namespace Tridiag

}	// namespace LNOT::BasicLinalg


#endif // LNOT_BASIC_LINALG_IMPL_HPP
