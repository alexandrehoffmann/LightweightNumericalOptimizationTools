#ifndef LNOT_MUTABLE_HPP
#define LNOT_MUTABLE_HPP

#include <LNOT/Fixed/Fixed.hpp>

namespace lnot
{

namespace detail
{

template <typename T>
struct MutableTraits
{
  using Type = T;
};

template <typename T, T value>
struct MutableTraits<Fixed<T, value>>
{
  using Type = T;
};

} // namespace detail

/**
 * @brief Mutable type metafunction
 * @tparam T The type for which to determine the mutable type
 * This alias template provides the mutable type corresponding to the given type T.
 * For Fixed types, it extracts the underlying type T.
 * For other types, it simply returns the type itself.
 */
template <typename T>
using Mutable = detail::MutableTraits<T>::Type;

} // namespace lnot

#endif // LNOT_MUTABLE_HPP
