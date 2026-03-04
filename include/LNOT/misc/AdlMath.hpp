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

template<typename T> inline constexpr std::remove_reference_t<T> pow(T&& x, T&& y) { using std::pow; return pow(std::forward<T>(x), std::forward<T>(y)); }

template<typename T> inline constexpr std::remove_reference_t<T> fma(T&& x, T&& y, T&& z) { using std::fma; return fma(std::forward<T>(x), std::forward<T>(y), std::forward<T>(z)); }

//

template<typename T> inline constexpr T fma(const T& x, const T& y, const T& z) { using std::fma; return fma(x, y, z); }
	
} // namespace LNOT::AdlMath

#endif // LNOT_ADL_MATH_HPP
