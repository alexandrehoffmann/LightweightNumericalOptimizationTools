#ifndef LNOT_TRUNCATED_CONJUGATE_GRADIENT_HPP
#define LNOT_TRUNCATED_CONJUGATE_GRADIENT_HPP

#include <LNOT/TRSSolvers/TRSSolverBase.hpp>

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
	
	void clearWorkSpace();
	
	void resizeWorkSpace(const Size newSize);
	
	template<typename HesOp, typename PrecOp, typename ASize> 
	Scalar solveImpl(const HesOp& H, const PrecOp& invB, const Scalar* g, const ASize size, const Scalar& delta, Scalar* x) requires (AreHessianOps<HesOp,PrecOp>::value and IsSize<ASize>::value);
	
	Scalar getErrorImpl        () const { return std::sqrt(m_precSqNormR); }
	Scalar getSquaredErrorImpl () const { return m_precSqNormR;            }
protected:
	LNOT_TRS_SOLVER_ATTRIBUTE
private:
	static Scalar getPolyMaxRoot(const Scalar a, const Scalar b, const Scalar c);

	Scalar m_precSqNormR = Scalar{};

	Scalar* m_z  = nullptr;
	Scalar* m_r  = nullptr;
	Scalar* m_p  = nullptr;
	Scalar* m_Hp = nullptr;
};

} // namespace LNOT

#include <LNOT/TRSSolvers/TruncatedConjugateGradient_impl.hpp>

#endif // LNOT_TRUNCATED_CONJUGATE_GRADIENT_HPP
