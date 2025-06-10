#ifndef LNOT_BASIC_LINALG_IMPL_HPP
#define LNOT_BASIC_LINALG_IMPL_HPP

#include <LNOT/BasicLinalg/BasicLinalg.hpp>

namespace LNOT
{
namespace BasicLinalg
{

template<typename Scalar, typename Size>
Scalar squaredNorm(const Scalar* __restrict__ x, const Size N) 
{ 
	Scalar res(0); 
	#pragma omp simd reduction(+:res)
	for (Size i=0; i!=N; ++i) { res += x[i]*x[i]; } 
	return res; 
}

#ifdef LNOT_WITH_BLAS
template<typename Size> float  squaredNorm(const float*  x, const Size N) { return cblas_sdot(blasint(N), x, 1, x, 1); }
template<typename Size> double squaredNorm(const double* x, const Size N) { return cblas_ddot(blasint(N), x, 1, x, 1); }
#endif // LNOT_WITH_BLAS 

template<typename Scalar, typename Size>
Scalar weightedSquaredNorm(const Scalar* __restrict__ x, const Scalar* __restrict__ w, const Size N)
{
	Scalar res(0);
	
	#pragma omp simd reduction(+:res)
	for (Size i=0; i!=N; ++i) {  res += x[i]*x[i]*w[i]; }
	
	return res;
}

template<typename Scalar, typename Size>
Scalar inner(const Scalar* __restrict__ x, const Scalar* __restrict__ y, const Size N)
{ 
	Scalar res(0); 
	#pragma omp simd reduction(+:res)
	for (Size i=0; i!=N; ++i) { res += x[i]*y[i]; } 
	return res; 
}

#ifdef LNOT_WITH_BLAS
template<typename Size> float  inner(const float*  x, const float*  y, const Size N) { return cblas_sdot(blasint(N), x, 1, y, 1); }
template<typename Size> double inner(const double* x, const double* y, const Size N) { return cblas_ddot(blasint(N), x, 1, y, 1); }
#endif // LNOT_WITH_BLAS 

template<typename Scalar, typename Size>
Scalar weightedInner(const Scalar* __restrict__ x, const Scalar* __restrict__ y, const Scalar* __restrict__ w, const Size N)
{ 
	Scalar res(0); 
	#pragma omp simd reduction(+:res)
	for (Size i=0; i!=N; ++i) { res += x[i]*y[i]*w[i]; } 
	return res; 
}

template<typename Scalar, typename Size>
void axpy(const Scalar alpha, const Scalar* __restrict__ x, const Size N, Scalar* __restrict__ y)
{
	#pragma omp simd 
	for (Size i=0; i!=N; ++i)
	{
		y[i] += alpha*x[i];
	}
}

#ifdef LNOT_WITH_BLAS
template<typename Size> void axpy(const float&  alpha, const float*  x, const Size N, float*  y) { cblas_saxpy(blasint(N), alpha, x, 1, y, 1); }
template<typename Size> void axpy(const double& alpha, const double* x, const Size N, double* y) { cblas_daxpy(blasint(N), alpha, x, 1, y, 1); }
#endif // LNOT_WITH_BLAS 

template<typename Scalar, typename Size>
void scal(const Scalar alpha, const Size N, Scalar* __restrict__ x)
{
	#pragma omp simd 
	for (Size i=0; i!=N; ++i)
	{
		x[i] *= alpha;
	}
}

#ifdef LNOT_WITH_BLAS
template<typename Size> void scal(const float  alpha, const Size N, float*  x) { cblas_sscal(blasint(N), alpha, x, 1); }
template<typename Size> void scal(const double alpha, const Size N, double* x) { cblas_dscal(blasint(N), alpha, x, 1); }
#endif // LNOT_WITH_BLAS 

template<typename Scalar, typename Size>
void symRk1Update(StorageOrder layout, UpLo uplo, const Scalar alpha, const Scalar* __restrict__ x, const Size N, Scalar* __restrict__ A)
{
	const Size iStride = layout == StorageOrder::ROW_MAJOR and uplo == UpLo::LOWER ? N : 1;
	const Size jStride = layout == StorageOrder::ROW_MAJOR and uplo == UpLo::LOWER ? 1 : N; 

	for (Size i=0; i!=N; ++i)
	{
		for (Size j=0; j!=(i+1); ++j)
		{
			A[i*iStride + j*jStride] += alpha*x[i]*x[j];
		}
	}
}

#ifdef LNOT_WITH_BLAS
template<typename Size> void symRk1Update(StorageOrder layout, UpLo uplo, const float  alpha, const float*  x, const Size N, float*  A) { cblas_ssyr(CBLAS_ORDER(layout), CBLAS_UPLO(uplo), blasint(N), alpha, x, 1, A, blasint(N)); }
template<typename Size> void symRk1Update(StorageOrder layout, UpLo uplo, const double alpha, const double* x, const Size N, double* A) { cblas_dsyr(CBLAS_ORDER(layout), CBLAS_UPLO(uplo), blasint(N), alpha, x, 1, A, blasint(N)); }
#endif // LNOT_WITH_BLAS 

namespace Tridiag
{

template<typename Scalar, typename Size>
Scalar norm1(const Scalar* __restrict__ alpha, const Scalar* __restrict__ beta, const Size N)
{
	Scalar res(0);
	#pragma omp simd reduction(+:res)
	for (Size i=0; i!=Size(N-1); ++i)
	{
		res += std::abs(alpha[i]) + 2*std::abs(beta[i]);
	}
	return res + std::abs(alpha[N-1]);
}

namespace LDLt
{

template<typename Scalar, typename Size>
bool compute(const Scalar* __restrict__ alpha, const Scalar* __restrict__ beta, const Size size, const Scalar shift, Scalar* __restrict__ invDelta, Scalar* __restrict__ l)
{
	if ((alpha[0] + shift) < std::numeric_limits<Scalar>::epsilon()) { return false; }
	
	invDelta[0] = Scalar(1) / (alpha[0] + shift);
	if (size == 1) { return true; } 
	l[0] = beta[0]*invDelta[0];
	
	for (Size i=1; i!=Size(size-1);++i)
	{
		const Scalar delta_i = alpha[i] + shift - beta[i-1]*l[i-1];
		if (delta_i < std::numeric_limits<Scalar>::epsilon()) { return false; }
		invDelta[i] = Scalar(1) / delta_i;
		l[i] = beta[i]*invDelta[i];
	}
	{ 
		const Size i = Size(size-1);
		const Scalar delta_i = alpha[i] + shift - beta[i-1]*l[i-1];
		if (delta_i < std::numeric_limits<Scalar>::epsilon()) { return false; }
		invDelta[i] = Scalar(1) / delta_i;
	}
	
	return true;
}

template<typename Scalar, typename Size>
void solve_L_e1(const Scalar* __restrict__ l, const Size size, const Scalar b1, Scalar* __restrict__ x)
{
	x[0] = b1;
	for (Size i=1; i!=size; ++i) { x[i] = -l[i-1]*x[i-1]; }
}

template<typename Scalar, typename Size>
void solve_L_inplace(const Scalar* __restrict__ l, const Size size, Scalar* __restrict__ x)
{
	for (Size i=1; i!=size; ++i) { x[i] -= l[i-1]*x[i-1]; }
}

template<typename Scalar, typename Size>
void solve_Lt_inplace(const Scalar* __restrict__ l, const Size size, Scalar* __restrict__ x)
{
	for (Size i=Size(size-2); i!=Size(-1); --i) {	x[i] -= l[i]*x[i+1]; }
}

template<typename Scalar, typename Size>
void solve_e1(const Scalar* __restrict__ invD, const Scalar* __restrict__ l, const Size size, const Scalar b1, Scalar* __restrict__ x)
{
	solve_L_e1(l, size, b1, x);     // Solve Lz = b
	for (Size i=0; i!=size; ++i) { x[i] *= invD[i]; }
	solve_Lt_inplace(l, size, x);   // Solve L^T x = y
}

} // namespace LDLt

} // namespace Tridiag

}	// namespace BasicLinalg
}


#endif // LNOT_BASIC_LINALG_IMPL_HPP
