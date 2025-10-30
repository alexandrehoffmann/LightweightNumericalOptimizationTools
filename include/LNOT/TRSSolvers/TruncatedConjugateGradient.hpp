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
	using Base = TRSSolverBase< TruncatedConjugateGradient<T> >;
public:
	using Scalar = typename Base::Scalar;
	using Size   = typename Base::Size;
	using Info   = typename Base::Info;

	template<typename HesOp, typename PrecOp> using AreHessianOps = typename Base::template AreHessianOps<HesOp,PrecOp>;
	
	void clearWorkSpace();
	
	void resizeWorkSpace(const Size newSize);
	
	template<typename HesOp, typename PrecOp> 
	void solve_impl(const HesOp& H, const PrecOp& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) requires (AreHessianOps<HesOp,PrecOp>::value);
	
	Scalar getError        () const { return std::sqrt(m_precSqNormR); }
	Scalar getSquaredError () const { return m_precSqNormR;            }
private:
	static Scalar getPolyMaxRoot(const Scalar a, const Scalar b, const Scalar c);

	Scalar m_precSqNormR = 0;

	Scalar* m_z  = nullptr;
	Scalar* m_r  = nullptr;
	Scalar* m_p  = nullptr;
	Scalar* m_Hp = nullptr;
};

using TruncatedConjugateGradientF = TruncatedConjugateGradient<float>;
using TruncatedConjugateGradientD = TruncatedConjugateGradient<double>;

} // namespace LNOT

#include <LNOT/TRSSolvers/TruncatedConjugateGradient_impl.hpp>

#endif // LNOT_TRUNCATED_CONJUGATE_GRADIENT_HPP
