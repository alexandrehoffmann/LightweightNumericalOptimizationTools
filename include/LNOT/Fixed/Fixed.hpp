#ifndef LNOT_FIXED_HPP
#define LNOT_FIXED_HPP

#include <LNOT/Settings.hpp>

namespace lnot
{
/**
 * @brief Fixed value wrapper
 * @tparam T The underlying type of the fixed value
 * @tparam value The fixed value of type T
 * This template struct represents a fixed value of type T.
 * The value is specified as a template parameter and can be accessed
 * via the VALUE static member or by converting the Fixed instance to T.
 */
template <typename T, T value>
struct Fixed
{
  using Type = T;

  static constexpr Type VALUE = value;

  // NOLINTNEXTLINE(cppcoreguidelines-explicit-constructor,-warnings-as-errors) : We want this to be implicit
  constexpr operator Type() const noexcept { return VALUE; }
};

/**
 * @brief Fixed value instance
 * @tparam T The underlying type of the fixed value
 * @tparam value The fixed value of type T
 * This variable template provides a convenient way to create
 * instances of the Fixed struct with the specified type and value.
 */
template <typename T, T value>
constexpr Fixed<T, value> FIXED = {};

/**
 * @brief Fixed value instance for default size used by LNOT
 * @tparam size The fixed size value
 * This variable template provides a convenient way to create
 * instances of the Fixed struct with the specified size.
 */
template <settings::DefaultSize size>
constexpr Fixed<settings::DefaultSize, size> FIXED_SIZE = {};

} // namespace lnot

#endif // LNOT_FIXED_HPP
