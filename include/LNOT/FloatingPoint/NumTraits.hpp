#ifndef LNOT_NUM_TRAITS_HPP
#define LNOT_NUM_TRAITS_HPP

#include <concepts>
#include <limits>

namespace LNOT
{

template<typename T> struct NumTraits;

template<std::floating_point T>
struct NumTraits<T>
{
	static constexpr T epsilon  = std::numeric_limits<T>::epsilon();
	static constexpr T infinity = std::numeric_limits<T>::infinity();
	static constexpr T lowest   = std::numeric_limits<T>::lowest();
};

} // namespace LNOT

#endif // LNOT_NUM_TRAITS_HPP
