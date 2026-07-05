#ifndef LNOT_LANCZOS_TRS_SOLVER_HPP
#define LNOT_LANCZOS_TRS_SOLVER_HPP

#include <LNOT/TRSSolvers/LanczosTRSSolverBase.hpp>

namespace LNOT
{

template<typename T> class LanczosTRSSolver;

template<typename T> 
struct TRSSolverTraits< LanczosTRSSolver<T> > 
{ 
	using Scalar = T; 
	using Size = unsigned int;
	
	static constexpr bool storeLanczosVectors = false;
};

template<typename T> 
class LanczosTRSSolver : public LanczosTRSSolverBase< LanczosTRSSolver<T> >
{
	using Self = LanczosTRSSolver<T>;
	
	friend class LanczosTRSSolverBase<Self>;
public:
	LNOT_DEFINE_TRS_SOLVER
	
	void resizeWorkSpace(const Size newSize);
protected:
	LNOT_TRS_SOLVER_ATTRIBUTE	
private:
	inline void resetBvImpl() {}
	
	inline Scalar* getBvCurrImpl() { return m_Bv_curr.get(); }
	
	inline Scalar* getBvPrevImpl() { return m_Bv_prev.get(); }
	
	inline const Scalar* getBvCurrImpl() const { return m_Bv_curr.get(); }
	
	inline const Scalar* getBvPrevImpl() const { return m_Bv_prev.get(); }
	
	template<typename ASize>
	void addBvNextImpl(const Scalar scale, const Scalar* unscaledBvNext, const ASize size) requires(isSize<ASize>);
	
	template<typename ASize>
	inline void reOrthonormalizeImpl(const ASize /* size */, Scalar* /* unscaledBvNext */ ) requires(isSize<ASize>) {}

	std::unique_ptr<Scalar[]> m_Bv_curr;
	std::unique_ptr<Scalar[]> m_Bv_prev;
};

} // namespace LNOT

#include <LNOT/TRSSolvers/LanczosTRSSolver_impl.hpp>

#endif // LNOT_LANCZOS_TRS_SOLVER_HPP
