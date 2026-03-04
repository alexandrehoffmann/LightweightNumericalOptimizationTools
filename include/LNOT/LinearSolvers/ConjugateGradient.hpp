#ifndef LNOT_CONJUGATE_GRADIENT_HPP
#define LNOT_CONJUGATE_GRADIENT_HPP

#include <LNOT/LinearSolvers/LinearSolverBase.hpp>
#include <LNOT/misc/AdlMath.hpp>

#include <cmath>

namespace LNOT
{

template<typename T> class ConjugateGradient;

template<std::floating_point T> struct LinearSolverTraits< ConjugateGradient<T> > { using Scalar = T;  using Size = unsigned int; };

template<typename T> 
class ConjugateGradient : public LinearSolverBase< ConjugateGradient<T> >
{
	using Self = ConjugateGradient<T>;
public:
	LNOT_DEFINE_LINEAR_SOLVER
	
	void clearWorkSpace();
	
	void resizeWorkSpace(const Size newSize); ///<  @brief reallocate internal memory if `newSize` > `Base::m_workCapacity`.

	template<typename HesOp, typename PrecOp, typename ASize, typename Bool>
	void solveImpl(const HesOp& H, const PrecOp& invB, const Scalar* g, const ASize size, const Bool solveInPlace, Scalar* x) requires (isHessianOp<HesOp> and isHessianOp<PrecOp> and isSize<ASize>);

	Scalar getErrorImpl() const { return AdlMath::sqrt(m_precSqNormR); }
	
	Scalar getSquaredErrorImpl() const { return m_precSqNormR; } 
protected:
	LNOT_LINEAR_SOLVER_ATTRIBUTE
private:
	Scalar m_precSqNormR = Scalar{};

	Scalar* m_z  = nullptr;
	Scalar* m_r  = nullptr;
	Scalar* m_p  = nullptr;
	Scalar* m_Hp = nullptr;
};

} // namespace LNOT

#include <LNOT/LinearSolvers/ConjugateGradient_impl.hpp>

#endif // LNOT_CONJUGATE_GRADIENT_HPP
