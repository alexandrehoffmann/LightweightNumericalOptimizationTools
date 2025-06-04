#ifndef LNOT_LINE_SEARCH_BASE_HPP
#define LNOT_LINE_SEARCH_BASE_HPP

#include <LNOT/Oracles/OracleBase.hpp>

#include <cstdio>

namespace LNOT
{

template<class LineSearch> struct LineSearchTraits;

template<class Derived>
class LineSearchBase
{
	using DerivedTraits = LineSearchTraits<Derived>;
public:
	using Scalar = typename DerivedTraits::Scalar;
	using Size   = typename DerivedTraits::Size;
	
	enum class Info {SUCCESS, FAILURE};

	const Derived& derived_cast() const { return static_cast<const Derived&>(*this); }
	      Derived& derived_cast()       { return static_cast<      Derived&>(*this); }
	
	LineSearchBase(const Size maxIt = 200000) : m_maxIt(maxIt) {}
	~LineSearchBase() { clearWorkSpace(); }
	
	void clearWorkSpace() { derived_cast().clearWorkSpace(); }
	
	template<FirstOrderOracle_concept Oracle>
	Scalar solve(const Scalar* x, const Scalar& fx, const Scalar* gx, const Scalar* s, Oracle& oracle) { derived_cast().solve(x, fx, gx, s, oracle); }
	
	Size getMaxIt   () const { return m_maxIt; }
	Size iterations () const { return m_nIt;   }
	Info info       () const { return m_info;  }
	
	void setMaxIt (const Size maxIt) { m_maxIt = maxIt; }
	
	void setOutput(std::FILE* out) { m_out = out; }
protected:
	Size m_maxIt;
	Size m_nIt;
	Info m_info;
	
	Size m_workCapacity = 0;
	
	std::FILE* m_out = nullptr;
};

template<class T> struct IsLineSearch : std::bool_constant< std::is_base_of<LineSearchBase<T>, T>::value > {};

template<class T> concept LineSearch_concept = IsLineSearch<T>::value;

} // 

#endif // LNOT_LINE_SEARCH_BASE_HPP
