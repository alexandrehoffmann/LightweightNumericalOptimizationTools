#ifndef LNOT_IS_MAYBE_FIXED_HPP
#define LNOT_IS_MAYBE_FIXED_HPP

#include <concepts>

#include <LNOT/Fixed/Fixed.hpp>
#include <LNOT/Fixed/Mutable.hpp>

namespace lnot
{

/**
 * @brief Concept to check if a type is a MaybeFixed value
 * @tparam MaybeFixedT The type to check if it is a MaybeFixed value
 * @tparam T The underlying type to compare against
 * This concept evaluates to true if the type, when passed through the Mutable metafunction,
 * results in the specified type T.
 */
template <typename MaybeFixedT, typename T>
concept Is_MaybeFixed = std::same_as<Mutable<MaybeFixedT>, T>;

// NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers) : Irrelevant
static_assert(Is_MaybeFixed<Fixed<int, 5>, int>, "Fixed<int, 5> should be recognized as MaybeFixed<int>.");
static_assert(Is_MaybeFixed<int, int>, "int should be recognized as MaybeFixed<int>.");
static_assert(!Is_MaybeFixed<Fixed<int, 5>, double>, "Fixed<int, 5> should not be recognized as MaybeFixed<double>.");
// NOLINTEND(cppcoreguidelines-avoid-magic-numbers)

} // namespace lnot

#endif // LNOT_IS_MAYBE_FIXED_HPP
