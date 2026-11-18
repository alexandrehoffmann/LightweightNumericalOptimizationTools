#ifndef LNOT_IS_FIXED_HPP
#define LNOT_IS_FIXED_HPP

#include <LNOT/Fixed/Fixed.hpp>

namespace lnot
{

namespace detail
{

template <typename T>
constexpr bool IS_FIXED = false;

template <typename T, T value>
constexpr bool IS_FIXED<Fixed<T, value>> = true;
} // namespace detail

/**
 * @brief Concept to check if a type is a Fixed value
 * @tparam T The type to check if it is a Fixed value
 * This concept evaluates to true if the type is an instance of the Fixed struct.
 */
template <typename T>
concept Is_Fixed = detail::IS_FIXED<T>;

// NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers) : False positive
static_assert(Is_Fixed<Fixed<int, 5>>, "Fixed<int, 5> should be recognized as a Fixed type.");
static_assert(!Is_Fixed<int>, "int should not be recognized as a Fixed type.");

} // namespace lnot

#endif // LNOT_IS_FIXED_HPP
