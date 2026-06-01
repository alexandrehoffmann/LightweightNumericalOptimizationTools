#ifndef LNOT_LANCZOS_SOLVER_HPP
#define LNOT_LANCZOS_SOLVER_HPP

#include <LNOT/LinearSolvers/LinearSolverBase.hpp>

#include <memory>

namespace LNOT
{

template<typename T> class LanczosSolver;

template<std::floating_point T> struct LinearSolverTraits< LanczosSolver<T> > { using Scalar = T;  using Size = unsigned int; };

template<typename T> 
class LanczosSolver : public LinearSolverBase< LanczosSolver<T> >
{
	using Self = LanczosSolver<T>;
public:
	LNOT_DEFINE_LINEAR_SOLVER
	
	void clearWorkSpace();
	
	void resizeWorkSpace(const Size newSize);  ///<  @brief reallocate internal memory if `newSize` > `Base::m_workCapacity`.
	
	template<typename HesOp, typename PrecOp, typename ASize, typename Bool>
	void solveImpl(const HesOp& H, const PrecOp& invB, const Scalar* g, const ASize size, const Bool solveInPlace, Scalar* x) requires(isHessianOp<HesOp> and isHessianOp<PrecOp> and isSize<ASize>);

	Scalar getErrorImpl() const { return m_precNormR;             }
	
	Scalar getSquaredErrorImpl() const { return m_precNormR*m_precNormR; }
protected:
	LNOT_LINEAR_SOLVER_ATTRIBUTE
private:
	Scalar m_precNormR = Scalar{};

	std::unique_ptr<Scalar[]> m_Bv_old;
	std::unique_ptr<Scalar[]> m_Bv;
	std::unique_ptr<Scalar[]> m_v;
	std::unique_ptr<Scalar[]> m_p;
	std::unique_ptr<Scalar[]> m_w;
};

} // namespace LNOT

#include <LNOT/LinearSolvers/LanczosSolver_impl.hpp>

#endif // LNOT_LANCZOS_SOLVER_HPP
