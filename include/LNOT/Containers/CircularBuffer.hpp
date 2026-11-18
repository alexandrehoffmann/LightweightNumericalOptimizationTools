#ifndef LNOT_CIRCULAR_BUFFER_HPP
#define LNOT_CIRCULAR_BUFFER_HPP

#include <iterator>
#include <memory>

#include <LNOT/Containers/ArrayHandle.hpp>
#include <LNOT/Fixed/Mutable.hpp>
#include <LNOT/Settings.hpp>

namespace lnot
{
// NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic) : This class needs pointer arithmetic

/**
 * CircularBuffer is a fixed-size circular buffer container.
 * It supports push operations and provides iterators for traversal.
 * The buffer overwrites the oldest elements when full.
 * @tparam ScalarT The type of elements stored in the buffer.
 * @tparam SizeT The type used for indexing and capacity.
 */
template <typename ScalarT, typename SizeT = settings::DefaultSize>
class CircularBuffer
{
public:
  using Scalar = ScalarT; ///<  @brief The type of elements stored in the buffer.
  using Size = SizeT;     ///<  @brief The type used for indexing and capacity.
  using DataHandle = HostArrayHandle<Scalar, Size>;

  class Sentinel; ///<  @brief Sentinel class for the circular buffer.

  template <bool isConst>
  class Iterator; ///<  @brief Iterator class for the circular buffer.

  /**
   * @brief Constructs a circular buffer with the specified capacity.
   * @param p_capacity The maximum number of elements the buffer can hold.
   */
  constexpr explicit CircularBuffer(Size p_capacity) noexcept;

  /**
   * @brief Pushes a value into the circular buffer.
   * @param p_value The value to be pushed.
   * @return The offset at which the value was inserted.
   */
  constexpr Mutable<Size> push(Scalar p_value) noexcept { return emplace(std::move(p_value)); }

  /**
   * @brief Clears the circular buffer, removing all elements.
   */
  constexpr void clear() noexcept;

  /**
   * @brief Constructs an element in-place in the circular buffer.
   * @tparam Args The types of the arguments to forward to the element's constructor.
   * @param args The arguments to forward to the element's constructor.
   * @return The offset at which the element was inserted.
   */
  template <typename... Args>
  constexpr Mutable<Size> emplace(Args&&... args) noexcept;

  /**
   * @brief Gets the offset of the first element in the circular buffer.
   * @return The offset of the first element.
   */
  constexpr Mutable<Size> getFrontOffset() const noexcept { return begin().getOffset(); }

  /**
   * @brief Gets the offset of the last element in the circular buffer.
   * @return The offset of the last element.
   */
  constexpr Mutable<Size> getBackOffset() const noexcept { return std::prev(std::ranges::next(begin(), end())).getOffset(); }

  /**
   * @brief Gets the offset of the bound (end) of the circular buffer.
   * @return The offset of the bound.
   */
  constexpr Mutable<Size> getBoundOffset() const noexcept;

  /**
   * @brief Gets the capacity of the circular buffer.
   * @return The capacity of the circular buffer.
   */
  constexpr Size getCapacity() const noexcept { return m_capacity; }

  /**
   * @brief Gets the current size of the circular buffer.
   * @return The number of elements currently stored in the circular buffer.
   */
  constexpr Mutable<Size> getSize() const noexcept;

  /**
   * @brief Returns an iterator to the beginning of the circular buffer.
   * @return An iterator to the beginning of the circular buffer.
   */
  constexpr Iterator<false> begin() noexcept { return { m_isFull ? m_bound : beginOfAllocated(), beginOfAllocated(), endOfAllocated() }; }

  /**
   * @brief Returns a const iterator to the beginning of the circular buffer.
   * @return A const iterator to the beginning of the circular buffer.
   */
  constexpr Iterator<true> begin() const noexcept;

  /**
   * @brief Returns a sentinel representing the end of the circular buffer.
   * @return A sentinel representing the end of the circular buffer.
   * @note The sentinel does not point to an actual element, but rather represents the number of steps to reach the end of the circular
   * buffer.
   */
  constexpr Sentinel end() const noexcept { return Sentinel(getSize()); } // cppcheck-suppress functionStatic

private:
  /**
   * @brief Returns a pointer to the beginning of the allocated memory.
   * @return A pointer to the beginning of the allocated memory.
   */
  constexpr Scalar* beginOfAllocated() const noexcept { return m_data.cPtr(); }

  /**
   * @brief Returns a pointer to the end of the allocated memory.
   * @return A pointer to the end of the allocated memory.
   */
  constexpr Scalar* endOfAllocated() const noexcept { return beginOfAllocated() + m_capacity; }

  DataHandle m_data;
  Scalar* m_bound;
  Size m_capacity;
  bool m_isFull = false;
};

// NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)

/**
 * @brief Sentinel class for the circular buffer.
 *
 * Represents the end of the circular buffer for iteration purposes.
 */
template <typename ScalarT, typename SizeT>
class CircularBuffer<ScalarT, SizeT>::Sentinel
{
public:
  constexpr Sentinel() noexcept = default;
  constexpr Sentinel(const Sentinel&) noexcept = default;
  constexpr Sentinel(Sentinel&&) noexcept = default;

  /**
   * @brief Constructs a sentinel with the given offset.
   * @param p_i The offset representing the number of steps to reach the end of the circular buffer.
   */
  constexpr explicit Sentinel(Mutable<Size> p_i) noexcept;

  constexpr ~Sentinel() noexcept = default;

  // NOLINTNEXTLINE(cppcoreguidelines-explicit-constructor,-warnings-as-errors) : We want this to be implicit
  constexpr operator Mutable<Size>() const noexcept { return m_i; }

  constexpr Sentinel& operator=(const Sentinel&) noexcept = default;

  constexpr Sentinel& operator=(Sentinel&&) noexcept = default;

private:
  Mutable<Size> m_i;
};

template <typename ScalarT, typename SizeT>
template <bool isConst>
class CircularBuffer<ScalarT, SizeT>::Iterator
{
public:
  template <typename U>
  using MaybeConst = std::conditional_t<isConst, const U, U>;

  /// @brief The underlying pointer type for the iterator.
  using InnerIt = MaybeConst<Scalar>*;
  /// @brief The value type of the elements pointed to by the iterator.
  using value_type = std::iter_value_t<InnerIt>;
  /// @brief The reference type of the elements pointed to by the iterator.
  using reference = std::iter_reference_t<InnerIt>;
  /// @brief The difference type used for the iterator.
  using difference_type = Mutable<Size>;
  /// @brief The iterator category tag, mainly used for iterator traits.
  using iterator_tag = std::random_access_iterator_tag;

  constexpr Iterator() noexcept = default;
  constexpr Iterator(const Iterator&) noexcept = default;
  constexpr Iterator(Iterator&&) noexcept = default;

  /**
   * @brief Constructs an iterator with the given current, begin, and bound pointers.
   * @param p_curr The current position of the iterator.
   * @param p_begin A pointer to the beginning of the allocated memory for the circular buffer.
   * @param p_bound A pointer to the end of the allocated memory for the circular buffer.
   */
  constexpr Iterator(InnerIt p_curr, InnerIt p_begin, InnerIt p_bound) noexcept;

  constexpr ~Iterator() noexcept = default;

  constexpr Iterator& operator=(const Iterator&) noexcept = default;

  constexpr Iterator& operator=(Iterator&&) noexcept = default;

  /**
   * @brief Returns the current offset of the iterator relative to the beginning of the allocated memory for the circular buffer.
   * @return The offset of the current element.
   */
  constexpr difference_type getOffset() const noexcept { return static_cast<difference_type>(std::distance(m_begin, m_curr)); }

  constexpr reference operator*() const noexcept { return *m_curr; }

  constexpr reference operator[](difference_type p_diff) const noexcept { return *(m_curr + p_diff); }

  constexpr Iterator& operator++() noexcept;

  constexpr Iterator& operator--() noexcept;

  constexpr Iterator operator++(int) noexcept;

  constexpr Iterator operator--(int) noexcept;

  constexpr Iterator& operator+=(difference_type p_diff) noexcept { return advance(p_diff); }

  constexpr Iterator& operator-=(difference_type p_diff) noexcept { return advance(-p_diff); }

  friend constexpr Iterator operator+(Iterator p_it, difference_type p_diff) noexcept { return Iterator(std::move(p_it)).advance(p_diff); }

  friend constexpr Iterator operator-(Iterator p_it, difference_type p_diff) noexcept { return Iterator(std::move(p_it)).advance(-p_diff); }

  friend constexpr Iterator operator+(difference_type p_diff, Iterator p_it) noexcept { return Iterator(std::move(p_it)).advance(p_diff); }

  friend constexpr bool operator==(const Iterator& p_lhs, const Iterator& p_rhs) noexcept = default;

  friend constexpr bool operator==(const Iterator& p_it, Sentinel p_bound) noexcept { return p_it.m_i == p_bound; }

  friend constexpr auto operator<=>(const Iterator& p_lhs, const Iterator& p_rhs) noexcept { return p_lhs.m_i <=> p_rhs.m_i; }

  friend constexpr difference_type operator-(const Iterator& p_lhs, const Iterator& p_rhs) noexcept { return p_lhs.m_i - p_rhs.m_i; }

  friend constexpr difference_type operator-(const Iterator& p_lhs, Sentinel p_bound) noexcept { return p_lhs.m_i - p_bound; }

  friend constexpr difference_type operator-(Sentinel p_bound, const Iterator& p_rhs) noexcept { return p_bound - p_rhs.m_i; }

private:
  constexpr Iterator& advance(difference_type p_diff) noexcept;

  InnerIt m_curr = {};
  InnerIt m_begin = {};
  InnerIt m_bound = {};
  difference_type m_i = {};
};

extern template class CircularBuffer<settings::DefaultScalar, settings::DefaultSize>;

static_assert(std::sized_sentinel_for<CircularBuffer<settings::DefaultScalar, settings::DefaultSize>::Sentinel,
                                      CircularBuffer<settings::DefaultScalar, settings::DefaultSize>::Iterator<false>>);
static_assert(std::sized_sentinel_for<CircularBuffer<settings::DefaultScalar, settings::DefaultSize>::Sentinel,
                                      CircularBuffer<settings::DefaultScalar, settings::DefaultSize>::Iterator<true>>);
static_assert(std::random_access_iterator<CircularBuffer<settings::DefaultScalar, settings::DefaultSize>::Iterator<false>>);
static_assert(std::random_access_iterator<CircularBuffer<settings::DefaultScalar, settings::DefaultSize>::Iterator<true>>);

} // namespace lnot

#include <LNOT/Containers/CircularBuffer_impl.hpp>

#endif // LNOT_CIRCULAR_BUFFER_HPP
