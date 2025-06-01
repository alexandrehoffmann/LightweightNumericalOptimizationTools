#ifndef LNOT_BASIC_LINALG_HPP
#define LNOT_BASIC_LINALG_HPP

#include <cmath>

namespace LNOT
{	
namespace BasicLinalg
{

template<typename Scalar, typename Size> Scalar squaredNorm(const Scalar* __restrict__ x, const Size N);
template<typename Scalar, typename Size> Scalar        norm(const Scalar* __restrict__ x, const Size N) { return std::sqrt(squaredNorm(x, N)); }

template<typename Scalar, typename Size> Scalar weightedSquaredNorm(const Scalar* __restrict__ x, const Scalar* __restrict__ w, const Size N);
template<typename Scalar, typename Size> Scalar        weightedNorm(const Scalar* __restrict__ x, const Scalar* __restrict__ w, const Size N) { return std::sqrt(weightedSquaredNorm(x, w, N)); }

template<typename Scalar, typename Size> Scalar         inner(const Scalar* __restrict__ x, const Scalar* __restrict__ y, const Size N);
template<typename Scalar, typename Size> Scalar weightedInner(const Scalar* __restrict__ x, const Scalar* __restrict__ y, const Scalar* __restrict__ w, const Size N);

template<typename Scalar, typename Size> void axpy(const Scalar& __restrict__ alpha, const Scalar* __restrict__ x, const Size N, Scalar* __restrict__ y);
template<typename Scalar, typename Size> void scal(const Scalar& __restrict__ alpha, const Size N, Scalar* __restrict__ x);

template<typename Scalar, typename Size> void symRk1Update(const Scalar& __restrict__ alpha, const Scalar* __restrict__ x, const Size N, Scalar* __restrict__ A);

namespace Tidiag
{

template<typename Scalar, typename Size> Scalar norm1(const Scalar* __restrict__ alpha, const Scalar* __restrict__ beta, const Size size);

namespace LDLt
{

template<typename Scalar, typename Size> bool compute(const Scalar* __restrict__ alpha, const Scalar* __restrict__ beta, const Size size, const Scalar shift, Scalar* __restrict__ invDelta, Scalar* __restrict__ l);
template<typename Scalar, typename Size> void solve_e1(const Scalar* __restrict__ invD, const Scalar* __restrict__ l, const Size size, const Scalar b1, Scalar* __restrict__ x);
template<typename Scalar, typename Size> void solve_L_e1(const Scalar* __restrict__ l, const Size size, const Scalar& __restrict__ b1, Scalar* __restrict__ x);
template<typename Scalar, typename Size> void solve_L_inplace(const Scalar* __restrict__ l, const Size size, Scalar* __restrict__ x);
template<typename Scalar, typename Size> void solve_Lt_inplace(const Scalar* __restrict__ l, const Size size, Scalar* __restrict__ x);
}

} // namespace Tridiag
} // namespace BasicLinalg
}

#include <LNOT/BasicLinalg_impl.hpp>

#endif // LNOT_BASIC_LINALG_HPP
