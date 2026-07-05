#ifndef LNOT_LINE_SEARCH_BASE_HPP
#define LNOT_LINE_SEARCH_BASE_HPP

#include <LNOT/Oracles/OracleBase.hpp>

#include <cstdio>

namespace LNOT
{

template<class LineSearch> struct LineSearchTraits;

template<class Derived>
class LineSearchBase : public CRTPBase<Derived>
{
	using DerivedTraits = LineSearchTraits<Derived>;
public:	
	using Scalar = typename DerivedTraits::Scalar; ///<  @brief The scalar type used in computations (e.g., float, double)
	using Size   = typename DerivedTraits::Size;   ///<  @brief The size type used for indexing and loop counters
	
	enum class Info {SUCCESS, FAILURE}; ///<  @brief Enumeration indicating solver termination status.
	
	LineSearchBase(const Size maxIt = 200000) : m_maxIt(maxIt) {}
	
	template<CFirstOrderOracle Oracle>
	Scalar solve(const Scalar* x, const Scalar& fx, const Scalar* gx, const Scalar* s, Oracle& oracle) { return this->derived().solveImpl(x, fx, gx, s, oracle); }
	
	Size getMaxIt   () const { return m_maxIt; }
	Size iterations () const { return m_nIt;   }
	Info info       () const { return m_info;  }
	
	void setMaxIt (const Size maxIt) { m_maxIt = maxIt; }
	
	void setOutput(std::FILE* out) { m_out = out; }
protected:
	Size m_maxIt;
	Size m_nIt;
	Info m_info;
	
	Size m_workCapacity = Size{};
	
	std::FILE* m_out = nullptr;
};

#define LNOT_DEFINE_LINESEARCH \
	using Base   = LineSearchBase<Self>; \
	using Size   = typename Base::Size; \
	using Scalar = typename Base::Scalar; \
	using Info   = typename Base::Info; \
	
#define LNOT_LINESEARCH_ATTRIBUTE \
	using Base::m_maxIt; \
	using Base::m_nIt; \
	using Base::m_info; \
	using Base::m_workCapacity; \
	using Base::m_out; \

template<class LineSearch> concept CLineSearch = std::derived_from<LineSearch, LineSearchBase<LineSearch>>;

} // 

#endif // LNOT_LINE_SEARCH_BASE_HPP
