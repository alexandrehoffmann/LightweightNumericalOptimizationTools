#ifndef LNOT_CONJUGATE_GRADIENT_HPP
#define LNOT_CONJUGATE_GRADIENT_HPP

#include <LNOT/LinearSolvers/LinearSolverBase.hpp>

#include <cmath>

namespace LNOT
{

template<typename T> class ConjugateGradient;

template<typename T> struct LinearSolverTraits< ConjugateGradient<T> > { using Scalar = T;  using Size = unsigned int; };

template<typename T> 
class ConjugateGradient : public LinearSolverBase< ConjugateGradient<T> >
{
	using Base = LinearSolverBase< ConjugateGradient<T> >;
public:
	using Scalar = typename Base::Scalar; ///<  @brief The scalar type used in computations (e.g., float, double)
	using Size   = typename Base::Size;   ///<  @brief The size type used for indexing and loop counters
	using Info   = typename Base::Info;   ///<  @brief Enumeration indicating solver termination status.

	template<typename HesOp, typename PrecOp> using AreHessianOps = typename Base::template AreHessianOps<HesOp,PrecOp>; ///<  @brief Trait to check if two types are both valid Hessian operators.
	
	void clearWorkSpace();
	
	void resizeWorkSpace(const Size newSize); ///<  @brief reallocate internal memory if `newSize` > `Base::m_workCapacity`.

	template<typename HesOp, typename PrecOp, bool solveInPlace> 
	void solve_impl(const HesOp& H, const PrecOp& invB, const Scalar* g, const Size size, std::bool_constant<solveInPlace>, Scalar* x) requires (AreHessianOps<HesOp,PrecOp>::value);

	Scalar getError        () const { return std::sqrt(m_precSqNormR); } 
	Scalar getSquaredError () const { return m_precSqNormR;            } 
private:
	Scalar m_precSqNormR = 0;

	Scalar* m_z  = nullptr;
	Scalar* m_r  = nullptr;
	Scalar* m_p  = nullptr;
	Scalar* m_Hp = nullptr;
};

using ConjugateGradientF = ConjugateGradient<float>;
using ConjugateGradientD = ConjugateGradient<double>;

} // namespace LNOT

#include <LNOT/LinearSolvers/ConjugateGradient_impl.hpp>

#endif // LNOT_CONJUGATE_GRADIENT_HPP
