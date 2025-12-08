#ifndef LNOT_LIMITED_MEMORY_SOLVER_BASE_HPP
#define LNOT_LIMITED_MEMORY_SOLVER_BASE_HPP

namespace LNOT
{

template<typename Size>
class LimitedMemorySolverBase
{
public:
	Size getMemory() const { return m_memory; }
	
	void setMemory(const Size memory) { m_memory = memory; }
protected:
	Size m_memory = Size(5);
};

} // namespace LNOT

#define LNOT_DEFINE_LIMITED_MEMORY_SOLVER \
	using LMBase = LimitedMemorySolverBase<Size>; \

#define LNOT_LIMITED_MEMORY_SOLVER_ATTRIBUTE \
	using LMBase::m_memory ;\

#endif // LNOT_LIMITED_MEMORY_SOLVER_BASE_HPP
