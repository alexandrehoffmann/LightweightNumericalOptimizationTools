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
	using Self = ConjugateGradient<T>;
public:
	LNOT_DEFINE_LINEAR_SOLVER
	
	void clearWorkSpace();
	
	void resizeWorkSpace(const Size newSize); ///<  @brief reallocate internal memory if `newSize` > `Base::m_workCapacity`.

	template<typename HesOp, typename PrecOp, typename ASize, typename Bool> requires(BIC::IsFixed<Bool>::value)
	void solveImpl(const HesOp& H, const PrecOp& invB, const Scalar* g, const ASize size, Bool solveInPlace, Scalar* x) requires (AreHessianOps<HesOp,PrecOp>::value and IsSize<ASize>::value);

	Scalar getErrorImpl        () const { return std::sqrt(m_precSqNormR); } 
	Scalar getSquaredErrorImpl () const { return m_precSqNormR;            } 
protected:
	LNOT_LINEAR_SOLVER_ATTRIBUTE
private:
	Scalar m_precSqNormR = 0;

	Scalar* m_z  = nullptr;
	Scalar* m_r  = nullptr;
	Scalar* m_p  = nullptr;
	Scalar* m_Hp = nullptr;
};

} // namespace LNOT

#include <LNOT/LinearSolvers/ConjugateGradient_impl.hpp>

#endif // LNOT_CONJUGATE_GRADIENT_HPP
