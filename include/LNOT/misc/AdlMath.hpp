#ifndef LNOT_ADL_MATH_HPP
#define LNOT_ADL_MATH_HPP

#include <cmath>
#include <utility>

namespace LNOT::AdlMath
{
	
template<typename T> inline constexpr std::remove_reference_t<T> abs(T&& x) { using std::abs; return abs(std::forward<T>(x)); }

template<typename T> inline constexpr std::remove_reference_t<T> sqrt(T&& x) { using std::sqrt; return sqrt(std::forward<T>(x)); }

template<typename T> inline constexpr std::remove_reference_t<T> pow(T&& x) { using std::pow; return pow(std::forward<T>(x)); }

template<typename T> inline constexpr std::remove_reference_t<T> floor(T&& x) { using std::floor; return floor(std::forward<T>(x)); }

template<typename T> inline constexpr std::remove_reference_t<T> log10(T&& x) { using std::log10; return log10(std::forward<T>(x)); }

template<typename T> inline constexpr bool isfinite(T&& x) { using std::isfinite; return isfinite(std::forward<T>(x)); }
	
}

#endif // LNOT_ADL_MATH_HPP
