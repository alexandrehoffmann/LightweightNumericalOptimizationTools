#ifndef LNOT_CIRCULAR_BUFFER_IMPL_HPP
#define LNOT_CIRCULAR_BUFFER_IMPL_HPP

#include <LNOT/Containers/CircularBuffer.hpp>

namespace lnot
{
// NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic) : This class needs pointer arithmetic

// CircularBuffer //

template <typename ScalarT, typename SizeT>
constexpr CircularBuffer<ScalarT, SizeT>::CircularBuffer(Size p_capacity) noexcept
    : m_data(DataHandle::allocate(p_capacity))
    , m_bound(m_data.cPtr())
    , m_capacity(std::move(p_capacity))
{
}

template <typename ScalarT, typename SizeT>
template <typename... Args>
constexpr auto CircularBuffer<ScalarT, SizeT>::emplace(Args&&... args) noexcept -> Mutable<Size>
{
  const Mutable<Size> pos = static_cast<Mutable<Size>>(std::distance(beginOfAllocated(), m_bound));

  std::destroy_at(m_bound);
  std::construct_at(m_bound, std::forward<Args>(args)...);
  ++m_bound;

  if (m_bound == endOfAllocated())
  {
    m_isFull = true;
    m_bound = beginOfAllocated();
  }

  return pos;
}

template <typename ScalarT, typename SizeT>
constexpr void CircularBuffer<ScalarT, SizeT>::clear() noexcept
{
  m_isFull = false;
  m_bound = beginOfAllocated();
}

template <typename ScalarT, typename SizeT>
constexpr auto CircularBuffer<ScalarT, SizeT>::getBoundOffset() const noexcept -> Mutable<Size>
{
  return static_cast<Mutable<Size>>(std::distance(beginOfAllocated(), m_bound));
}

template <typename ScalarT, typename SizeT>
constexpr auto CircularBuffer<ScalarT, SizeT>::getSize() const noexcept -> Mutable<Size>
{
  return m_isFull ? Mutable<Size>(m_capacity) : static_cast<Mutable<Size>>(std::distance(beginOfAllocated(), m_bound));
}

template <typename ScalarT, typename SizeT>
constexpr auto CircularBuffer<ScalarT, SizeT>::begin() const noexcept -> Iterator<true>
{
  return { m_isFull ? m_bound : beginOfAllocated(), beginOfAllocated(), endOfAllocated() };
}

// CircularBuffer::Sentinel //

template <typename ScalarT, typename SizeT>
constexpr CircularBuffer<ScalarT, SizeT>::Sentinel::Sentinel(Mutable<Size> p_i) noexcept
    : m_i(std::move(p_i))
{
}

// CircularBuffer::Iterator //

template <typename ScalarT, typename SizeT>
template <bool isConst>
constexpr CircularBuffer<ScalarT, SizeT>::Iterator<isConst>::Iterator(InnerIt p_curr, InnerIt p_begin, InnerIt p_bound) noexcept
    : m_curr(std::move(p_curr))
    , m_begin(std::move(p_begin))
    , m_bound(std::move(p_bound))
{
}

template <typename ScalarT, typename SizeT>
template <bool isConst>
constexpr auto CircularBuffer<ScalarT, SizeT>::Iterator<isConst>::operator++() noexcept -> Iterator&
{
  ++m_curr;
  ++m_i;
  if (m_curr == m_bound)
  {
    m_curr = m_begin;
  }
  return *this;
}

template <typename ScalarT, typename SizeT>
template <bool isConst>
constexpr auto CircularBuffer<ScalarT, SizeT>::Iterator<isConst>::operator--() noexcept -> Iterator&
{
  if (m_curr == m_begin)
  {
    m_curr = m_bound;
  }
  --m_curr;
  --m_i;

  return *this;
}

template <typename ScalarT, typename SizeT>
template <bool isConst>
constexpr auto CircularBuffer<ScalarT, SizeT>::Iterator<isConst>::operator++(int) noexcept -> Iterator
{
  Iterator tmp(*this);
  ++(*this);
  return tmp;
}

template <typename ScalarT, typename SizeT>
template <bool isConst>
constexpr auto CircularBuffer<ScalarT, SizeT>::Iterator<isConst>::operator--(int) noexcept -> Iterator
{
  Iterator tmp(*this);
  --(*this);
  return tmp;
}

template <typename ScalarT, typename SizeT>
template <bool isConst>
constexpr auto CircularBuffer<ScalarT, SizeT>::Iterator<isConst>::advance(difference_type p_diff) noexcept -> Iterator&
{
  if (p_diff >= 0)
  {
    difference_type remaining = static_cast<difference_type>(std::ranges::advance(m_curr, p_diff, m_bound));
    while (remaining > 0)
    {
      m_curr = m_begin;
      remaining = static_cast<difference_type>(std::ranges::advance(m_curr, remaining, m_bound));
    }
    if (m_curr == m_bound)
    {
      m_curr = m_begin;
    }
    m_i += p_diff;
  }
  else
  {
    if (m_curr == m_begin and p_diff != 0)
    {
      m_curr = m_bound;
    }
    difference_type remaining = static_cast<difference_type>(std::ranges::advance(m_curr, p_diff, m_begin));
    while (remaining < 0)
    {
      m_curr = m_bound;
      remaining = static_cast<difference_type>(std::ranges::advance(m_curr, remaining, m_begin));
    }
    m_i += p_diff;
  }
  return *this;
}

// NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)

} // namespace lnot

#endif // LNOT_CIRCULAR_BUFFER_IMPL_HPP
