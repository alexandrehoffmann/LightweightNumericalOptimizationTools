#ifndef LNOT_DIAGONAL_PRECONDITIONER_OP_IMPL_HPP
#define LNOT_DIAGONAL_PRECONDITIONER_OP_IMPL_HPP

#include <LNOT/BasicLinalg/DiagonalPreconditionerOp.hpp>

namespace LNOT
{

//// explicit template instanciations ////

extern template class DiagonalPreconditionerOp<float>;
extern template class DiagonalPreconditionerOp<double>;
extern template class DiagonalPreconditionerOp<long double>;

//// method implementations ////

template<std::floating_point Scalar>
void DiagonalPreconditionerOp<Scalar>::recompute(const Scalar* A) 
	: m_invD(new Scalar[size])
	, m_size(size) 
{ 
	#pragma omp simd
	for (Size i=0; i!=m_size; ++i) { m_invD[i] = Scalar(1) / A[i*(size + 1)]; } 
}

template<std::floating_point Scalar>
void DiagonalPreconditionerOp<Scalar>::operator() (const Scalar* x, Scalar* invDx) const
{
	#pragma omp simd
	for (Size i=0; i!=m_size; ++i) { invDx[i] = m_invD[i]*x[i]; }
}

} // namespace LNOT


#endif // LNOT_DIAGONAL_PRECONDITIONER_OP_IMPL_HPP

