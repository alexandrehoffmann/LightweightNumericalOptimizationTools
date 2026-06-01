#ifndef LNOT_LANCZOS_FULL_ORTH_TRS_SOLVER_HPP
#define LNOT_LANCZOS_FULL_ORTH_TRS_SOLVER_HPP

#include <LNOT/TRSSolvers/LanczosTRSSolverBase.hpp>

#include <vector>

namespace LNOT
{

template<typename T, size_t nPass> class LanczosFullOrthTRSSolver;

template<typename T, size_t nPass> 
struct TRSSolverTraits< LanczosFullOrthTRSSolver<T, nPass> > 
{ 
	using Scalar = T; 
	using Size = unsigned int;
	
	static constexpr bool storeLanczosVectors = true;
};

template<typename T, size_t nPass = 1> 
class LanczosFullOrthTRSSolver : public LanczosTRSSolverBase< LanczosFullOrthTRSSolver<T, nPass> >
{
	using Self = LanczosFullOrthTRSSolver<T, nPass>;
	
	friend class LanczosTRSSolverBase<Self>;
public:
	LNOT_DEFINE_TRS_SOLVER
	
	void resizeWorkSpace(const Size newSize);
protected:
	LNOT_TRS_SOLVER_ATTRIBUTE	
private:
	inline void resetBvImpl() { m_prev = {}; }
	
	inline Scalar* getBvCurrImpl() { return getBvImpl(m_prev + 1); }
	
	inline Scalar* getBvPrevImpl() { return getBvImpl(m_prev); }
	
	inline Scalar* getBvImpl(const Size i) { return m_BV.data() + i*m_vecSize; }
	
	inline const Scalar* getBvCurrImpl() const { return getBvImpl(m_prev + 1); }
	
	inline const Scalar* getBvPrevImpl() const { return getBvImpl(m_prev); }
	
	inline const Scalar* getBvImpl(const Size i) const { return m_BV.data() + i*m_vecSize; }
	
	template<typename ASize>
	void addBvNextImpl(const Scalar scale, const Scalar* unscaledBvNext, const ASize size) requires(IsSize<ASize>::value);
	
	template<typename ASize>
	void reOrthonormalizeImpl(const ASize size, Scalar* unscaledBvNext) requires(IsSize<ASize>::value);

	std::vector<Scalar> m_BV;
	Size                m_prev = {};
	Size                m_vecSize;
};

} // namespace LNOT

#include <LNOT/TRSSolvers/LanczosFullOrthTRSSolver_impl.hpp>

#endif // LNOT_LANCZOS_FULL_ORTH_TRS_SOLVER_HPP
