#ifndef LNOT_LANCZOS_SOLVER_HPP
#define LNOT_LANCZOS_SOLVER_HPP

#include <LNOT/LinearSolvers/LinearSolverBase.hpp>

namespace LNOT
{

template<typename T> class LanczosSolver;

template<typename T> struct LinearSolverTraits< LanczosSolver<T> > { using Scalar = T;  using Size = unsigned int; };

template<typename T> 
class LanczosSolver : public LinearSolverBase< LanczosSolver<T> >
{
	using Self = LanczosSolver<T>;
public:
	LNOT_DEFINE_LINEAR_SOLVER
	
	void clearWorkSpace();
	
	void resizeWorkSpace(const Size newSize);  ///<  @brief reallocate internal memory if `newSize` > `Base::m_workCapacity`.
	
	template<typename HesOp, typename PrecOp, typename ASize, typename Bool> requires(BIC::IsFixed<Bool>::value)
	void solveImpl(const HesOp& H, const PrecOp& invB, const Scalar* g, const ASize size, const Bool solveInPlace, Scalar* x) requires (AreHessianOps<HesOp,PrecOp>::value and IsSize<ASize>::value);

	Scalar getErrorImpl        () const { return m_precNormR;             }
	Scalar getSquaredErrorImpl () const { return m_precNormR*m_precNormR; }
protected:
	LNOT_LINEAR_SOLVER_ATTRIBUTE
private:
	Scalar m_precNormR = Scalar{};

	Scalar* m_Bv_old = nullptr;
	Scalar* m_Bv     = nullptr;
	Scalar* m_v      = nullptr;
	Scalar* m_p      = nullptr;
	Scalar* m_w      = nullptr;
};

} // namespace LNOT

#include <LNOT/LinearSolvers/LanczosSolver_impl.hpp>

#endif // LNOT_LANCZOS_SOLVER_HPP
