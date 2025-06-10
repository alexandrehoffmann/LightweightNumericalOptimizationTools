#ifndef LNOT_BASIC_LINALG_HPP
#define LNOT_BASIC_LINALG_HPP

#include <cmath>

#ifdef LNOT_WITH_BLAS
#include <cblas.h>
#endif // LNOT_WITH_BLAS

namespace LNOT
{
	
#ifdef LNOT_WITH_BLAS
enum class StorageOrder { ROW_MAJOR = CblasRowMajor, COL_MAJOR = CblasColMajor}; ///<  @brief Enum for matrix storage order.
enum class UpLo         { LOWER = CblasLower, UPPER = CblasUpper};               ///<  @brief Enum that specifies whether the upper or lower triangular part of 
#else
enum class StorageOrder { ROW_MAJOR, COL_MAJOR }; ///<  @brief Enum for matrix storage order.
enum class UpLo         { LOWER, UPPER };         ///<  @brief Enum that specifies whether the upper or lower triangular part of 
#endif // LNOT_WITH_BLAS

namespace BasicLinalg
{

///  @brief Computes \f$\sum_{i=1}^N x_i^2\f$
template<typename Scalar, typename Size> Scalar squaredNorm(const Scalar* __restrict__ x, const Size N);
///  @brief Computes \f$\|x\|_{\ell^2}\f$
template<typename Scalar, typename Size> Scalar norm(const Scalar* __restrict__ x, const Size N) { return std::sqrt(squaredNorm(x, N)); } 
///  @brief Computes \f$\sum_{i=1}^N w_i x_i^2\f$
template<typename Scalar, typename Size> Scalar weightedSquaredNorm(const Scalar* __restrict__ x, const Scalar* __restrict__ w, const Size N); 
///  @brief Computes \f$\sqrt{\sum_{i=1}^N w_i x_i^2}\f$
template<typename Scalar, typename Size> Scalar weightedNorm(const Scalar* __restrict__ x, const Scalar* __restrict__ w, const Size N) { return std::sqrt(weightedSquaredNorm(x, w, N)); } 
///  @brief Computes \f$\sqrt{\sum_{i=1}^N x_i y_i}\f$
template<typename Scalar, typename Size> Scalar inner(const Scalar* __restrict__ x, const Scalar* __restrict__ y, const Size N);
///  @brief Computes \f$\sqrt{\sum_{i=1}^N x_i w_i y_i}\f$
template<typename Scalar, typename Size> Scalar weightedInner(const Scalar* __restrict__ x, const Scalar* __restrict__ y, const Scalar* __restrict__ w, const Size N);
///  @brief Performs \f$ y = y + \alpha x \f$
template<typename Scalar, typename Size> void axpy(const Scalar alpha, const Scalar* __restrict__ x, const Size N, Scalar* __restrict__ y);
///  @brief Performs \f$ x = \alpha x \f$
template<typename Scalar, typename Size> void scal(const Scalar alpha, const Size N, Scalar* __restrict__ x);
/**
 * @brief Performs a Symmetric Rank 1 update of \f$A = A + \alpha xx^T\f$
 * @param layout Specifies whether two-dimensional array storage is row-major or column-major
 * @param uplo Specifies whether the upper or lower triangular part of the array A is used
 * @param alpha Scaling factor
 * @param x Input vector
 * @param N Vector/Matrix size
 * @param A An array of size `N*N` that stores either the lower or upper triangular part of the matrix symmetric \f$A\f$.
 */
template<typename Scalar, typename Size> void symRk1Update(StorageOrder layout, UpLo uplo, const Scalar alpha, const Scalar* __restrict__ x, const Size N, Scalar* __restrict__ A);

namespace Tridiag
{
/**
 * @brief Compute the \f$L^1\f$ norm of a tridiagonal matrix \f$T\in\mathbb{R}^N\f$
 * @param alpha An array of size N  that represents the diagonal of \f$T\f$
 * @param beta An array of size N-1 that represents the sub-diagonal of \f$T\f$
 * @param N The size of the matrix
 */
template<typename Scalar, typename Size> Scalar norm1(const Scalar* __restrict__ alpha, const Scalar* __restrict__ beta, const Size N);

namespace LDLt
{
/**
 * @brief Perform \f$LDL^T\f$ factorization of a symmetric tridiagonal matrix with shift.
 * @param alpha An array of size N  that represents the diagonal of the matrix to be factorized
 * @param beta An array of size N-1  that represents the sub-diagonal of the matrix to be factorized
 * @param size Matrix size
 * @param shift Diagonal shift to apply
 * @param invDelta Output array of size N. Stores the inverse of D
 * @param l Output array of size N. Stores lower bidiagonal elements (unit diagonal assumed)
 * @return true if factorization succeeded
 */
template<typename Scalar, typename Size> bool compute(const Scalar* __restrict__ alpha, const Scalar* __restrict__ beta, const Size size, const Scalar shift, Scalar* __restrict__ invDelta, Scalar* __restrict__ l);
/**
 * @brief Solve \f$LDL^Tx = b\f$ with \f$b = (b1, 0, ..., 0)^T\f$ using precomputed inverse diagonal.
 * @param invD Input array of size N. Stores the inverse of D
 * @param l Input array of size N. Stores lower bidiagonal elements (unit diagonal assumed)
 * @param size Matrix size
 * @param b1 First entry of the RHS vector
 * @param x Output array of size N. Stores the solution
 */
template<typename Scalar, typename Size> void solve_e1(const Scalar* __restrict__ invD, const Scalar* __restrict__ l, const Size size, const Scalar b1, Scalar* __restrict__ x);
/**
 * @brief Solve \f$Lx = b\f$ where \f$b = (b1, 0, ..., 0)^T\f$.
 * @param l Input array of size N. Stores lower bidiagonal elements (unit diagonal assumed)
 * @param size Matrix size
 * @param b1 First entry of the RHS vector
 * @param x Output array of size N. Stores the solution
 */
template<typename Scalar, typename Size> void solve_L_e1(const Scalar* __restrict__ l, const Size size, const Scalar b1, Scalar* __restrict__ x);
/**
 * @brief In-place solution of \f$Lx = b\f$ for general b (overwrites x).
 * @param l Input array of size N. Stores lower bidiagonal elements (unit diagonal assumed)
 * @param size Matrix size
 * @param x Array of size N. On entry, stores b. On exit stores the solution
 */
template<typename Scalar, typename Size> void solve_L_inplace(const Scalar* __restrict__ l, const Size size, Scalar* __restrict__ x);
/**
 * @brief In-place solution of \f$L^Tx = b\f$ for general b (overwrites x).
 * @param l Input array of size N. Stores lower bidiagonal elements (unit diagonal assumed)
 * @param size Matrix size
 * @param x Array of size N. On entry, stores b. On exit stores the solution
 */
template<typename Scalar, typename Size> void solve_Lt_inplace(const Scalar* __restrict__ l, const Size size, Scalar* __restrict__ x);
}

} // namespace Tridiag
} // namespace BasicLinalg
}

#include <LNOT/BasicLinalg/BasicLinalg_impl.hpp>

#endif // LNOT_BASIC_LINALG_HPP
