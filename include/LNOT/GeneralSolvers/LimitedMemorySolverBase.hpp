#ifndef LNOT_LIMITED_MEMORY_SOLVER_BASE_HPP
#define LNOT_LIMITED_MEMORY_SOLVER_BASE_HPP

template<typename Size>
class LimitedMemorySolverBase
{
public:
	Size getMemory() const { return m_memory; }
	
	void setMemory(const Size memory) { m_memory = memory; }
protected:
	Size m_memory = Size(5);
};

#endif // LNOT_LIMITED_MEMORY_SOLVER_BASE_HPP
