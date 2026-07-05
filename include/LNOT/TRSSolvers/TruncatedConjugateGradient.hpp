#ifndef LNOT_TRUNCATED_CONJUGATE_GRADIENT_HPP
#define LNOT_TRUNCATED_CONJUGATE_GRADIENT_HPP

#include <LNOT/TRSSolvers/TRSSolverBase.hpp>
#include <LNOT/misc/AdlMath.hpp>

#include <memory>

namespace LNOT
{

template<typename T> class TruncatedConjugateGradient;

template<typename T> struct TRSSolverTraits< TruncatedConjugateGradient<T> > { using Scalar = T; using Size = unsigned int; };

template<typename T>
class TruncatedConjugateGradient : public TRSSolverBase< TruncatedConjugateGradient<T> >
{
	using Self = TruncatedConjugateGradient<T>;
public:
	LNOT_DEFINE_TRS_SOLVER
	
	void resizeWorkSpace(const Size newSize);
	
	template<typename HesOp, typename PrecOp, typename ASize> 
	Scalar solveImpl(const HesOp& H, const PrecOp& invB, const Scalar* g, const ASize size, const Scalar& delta, Scalar* x) requires (areHessianOps<HesOp,PrecOp> and isSize<ASize>);
	
	Scalar getErrorImpl() const { return AdlMath::sqrt(m_precSqNormR); }
	
	Scalar getSquaredErrorImpl () const { return m_precSqNormR; }
protected:
	LNOT_TRS_SOLVER_ATTRIBUTE
private:
	static Scalar getPolyMaxRoot(const Scalar a, const Scalar b, const Scalar c);

	Scalar m_precSqNormR = Scalar{};

	std::unique_ptr<Scalar[]> m_z;
	std::unique_ptr<Scalar[]> m_r;
	std::unique_ptr<Scalar[]> m_p;
	std::unique_ptr<Scalar[]> m_Hp;
};

} // namespace LNOT

#include <LNOT/TRSSolvers/TruncatedConjugateGradient_impl.hpp>

#endif // LNOT_TRUNCATED_CONJUGATE_GRADIENT_HPP
