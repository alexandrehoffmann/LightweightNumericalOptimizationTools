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
	using Base = LinearSolverBase< LanczosSolver<T> >;
public:
	using Scalar = typename Base::Scalar; ///<  @brief The scalar type used in computations (e.g., float, double)
	using Size   = typename Base::Size;   ///<  @brief The size type used for indexing and loop counters
	using Info   = typename Base::Info;   ///<  @brief Enumeration indicating solver termination status.

	template<typename HesOp, typename PrecOp> using AreHessianOps = typename Base::template AreHessianOps<HesOp,PrecOp>; ///<  @brief Trait to check if two types are both valid Hessian operators.
	template<typename ASize>                  using IsSize        = typename Base::template IsSize<ASize>;               ///<  @brief Trait to check if a type is either a `Size` or a `BIC::Fixed<Size, VALUE>`
	
	void clearWorkSpace();
	
	void resizeWorkSpace(const Size newSize);  ///<  @brief reallocate internal memory if `newSize` > `Base::m_workCapacity`.
	
	template<typename HesOp, typename PrecOp, typename ASize, bool solveInPlace> 
	void solveImpl(const HesOp& H, const PrecOp& invB, const Scalar* g, const ASize size, std::bool_constant<solveInPlace>, Scalar* x) requires (AreHessianOps<HesOp,PrecOp>::value and IsSize<ASize>::value);

	Scalar getErrorImpl        () const { return m_precNormR;             }
	Scalar getSquaredErrorImpl () const { return m_precNormR*m_precNormR; }
private:
	Scalar m_precNormR = 0;

	Scalar* m_Bv_old = nullptr;
	Scalar* m_Bv     = nullptr;
	Scalar* m_v      = nullptr;
	Scalar* m_p      = nullptr;
	Scalar* m_w      = nullptr;
};

using LanczosSolverF = LanczosSolver<float>;
using LanczosSolverD = LanczosSolver<double>;

} // namespace LNOT

#include <LNOT/LinearSolvers/LanczosSolver_impl.hpp>

#endif // LNOT_LANCZOS_SOLVER_HPP
