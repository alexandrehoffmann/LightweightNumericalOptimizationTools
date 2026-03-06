#ifndef LNOT_LANCZOS_TRS_SOLVER_BASE_HPP
#define LNOT_LANCZOS_TRS_SOLVER_BASE_HPP

#include <LNOT/TRSSolvers/TRSSolverBase.hpp>

#include <vector>

namespace LNOT
{

template<class Derived> class LanczosTRSSolverBase;

template<class Derived> 
struct TRSSolverTraits< LanczosTRSSolverBase<Derived> > 
{ 
	using Scalar = typename TRSSolverTraits<Derived>::Scalar; 
	using Size   = typename TRSSolverTraits<Derived>::Size; 
};

template<class Derived> 
class LanczosTRSSolverBase : public TRSSolverBase<Derived>
{
	using Self = Derived;
public:
	LNOT_DEFINE_TRS_SOLVER
	using CRTP = typename Base::CRTP;
	
	static constexpr bool storeLanczosVectors = TRSSolverTraits<Derived>::storeLanczosVectors;
	
	LanczosTRSSolverBase(const Size maxIt = 200000, const Scalar tol = NumTraits<Scalar>::epsilon, const Size maxItTr = 200000, const Scalar tolTr = AdlMath::sqrt(NumTraits<Scalar>::epsilon));
	
	template<typename HesOp, typename PrecOp, typename ASize> 
	Scalar solveImpl(const HesOp& H, const PrecOp& invB, const Scalar* g, const ASize size, const Scalar& delta, Scalar* x) requires(AreHessianOps<HesOp,PrecOp>::value and IsSize<ASize>::value);
	
	Scalar getErrorImpl        () const { return m_precNormR;             }
	Scalar getSquaredErrorImpl () const { return m_precNormR*m_precNormR; }
	
	Size getMaxItTr() const { return m_maxItTr; }
	
	void setMaxItTr(const Size maxItTr) { m_maxItTr = maxItTr; }
	
	Scalar getLambda() const { return m_lambda; } 
protected:
	LNOT_TRS_SOLVER_ATTRIBUTE
	
	void clearLanczosWorkSpace();
	
	void allocateLanczosWorkSpace();
private:
	bool solveBoundary(const Scalar& gamma, const Scalar& delta);
	
	inline void resetBv() { CRTP::derived().resetBvImpl(); }
	
	inline Scalar* getBvCurr() { return CRTP::derived().getBvCurrImpl(); }
	
	inline Scalar* getBvPrev() { return CRTP::derived().getBvPrevImpl(); }
	
	inline Scalar* getBv(const Size i) requires(storeLanczosVectors) { return CRTP::derived().getBvImpl(i); }
	
	inline const Scalar* getBvCurr() const { return CRTP::derived().getBvCurrImpl(); }
	
	inline const Scalar* getBvPrev() const { return CRTP::derived().getBvPrevImpl(); }
	
	inline Scalar* getBv(const Size i) const requires(storeLanczosVectors) { return CRTP::derived().getBvImpl(i); }
	
	template<typename ASize>
	inline void addBvNext(const Scalar scale, const Scalar* unscaledBvNext, const ASize size) requires(IsSize<ASize>::value) { CRTP::derived().addBvNextImpl(scale, unscaledBvNext, size); } 
	
	template<typename ASize>
	inline void reOrthonormalize(const ASize size, Scalar* unscaledBvNext) requires(IsSize<ASize>::value) { CRTP::derived().reOrthonormalizeImpl(size, unscaledBvNext); } 
	
	Size m_maxItTr;

	Scalar m_precNormR = Scalar{};
	Scalar m_lambda    = Scalar{};

	Scalar* m_v      = nullptr;
	Scalar* m_p      = nullptr;
	Scalar* m_Bp     = nullptr;
	Scalar* m_Hp     = nullptr;
	Scalar* m_w      = nullptr;
	
	std::vector<Scalar> m_alpha;
	std::vector<Scalar> m_beta;
	
	std::vector<Scalar> m_invD;
	std::vector<Scalar> m_l;
	std::vector<Scalar> m_h;
};
	
} // namespace LNOT

#include <LNOT/TRSSolvers/LanczosTRSSolverBase_impl.hpp>

#endif // LNOT_LANCZOS_TRS_SOLVER_BASE_HPP
