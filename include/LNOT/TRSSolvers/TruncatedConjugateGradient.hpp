#ifndef LNOT_TRUNCATED_CONJUGATE_GRADIENT_HPP
#define LNOT_TRUNCATED_CONJUGATE_GRADIENT_HPP

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

	template<typename HessianOp> using IsHessianOp = typename Base::template IsHessianOp<HessianOp>;
	
	void clearWorkSpace();
	
	template<typename Op> 
	void solve(const Op& H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) requires (IsHessianOp<Op>::value);
	
	Scalar getError        () const { return std::sqrt(m_sqNormR); }
	Scalar getSquaredError () const { return m_sqNormR;            }
private:
	static Scalar getPolyMaxRoot(const Scalar a, const Scalar b, const Scalar c);

	Scalar m_sqNormR = 0;

	Scalar* m_r  = nullptr;
	Scalar* m_p  = nullptr;
	Scalar* m_Hp = nullptr;
};

} // namespace LNOT

#include <LNOT/TRSSolvers/TruncatedConjugateGradient_impl.hpp>

#endif // LNOT_TRUNCATED_CONJUGATE_GRADIENT_HPP
