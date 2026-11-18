#ifndef LNOT_ARRAY_HANDLE_IMPL_HPP
#define LNOT_ARRAY_HANDLE_IMPL_HPP

#include <LNOT/Containers/ArrayHandle.hpp>

namespace lnot
{

template <class AllocatorT>
constexpr ArrayHandle<AllocatorT>::ArrayHandle(Pointer p_handle) noexcept
    : m_handle(std::move(p_handle))
{
}

template <class AllocatorT>
constexpr ArrayHandle<AllocatorT>::ArrayHandle(ArrayHandle&& p_other) noexcept
    : m_handle(std::move(p_other.m_handle))
{
  p_other.m_handle = nullptr;
}

template <class AllocatorT>
constexpr ArrayHandle<AllocatorT>& ArrayHandle<AllocatorT>::operator=(ArrayHandle&& p_other) noexcept
{
  free();
  m_handle = std::move(p_other.m_handle);
  p_other.m_handle = nullptr;

  return *this;
}

template <class AllocatorT>
constexpr ArrayHandle<AllocatorT>& ArrayHandle<AllocatorT>::operator=(std::nullptr_t /* p_handle */) noexcept
{
  free();
  return *this;
}

template <class AllocatorT>
void ArrayHandle<AllocatorT>::free() noexcept
{
  if (m_handle != nullptr)
  {
    Allocator::free(m_handle);
    m_handle = nullptr;
  }
}

} // namespace lnot

#endif // LNOT_ARRAY_HANDLE_IMPL_HPP
