#ifndef LNOT_ARRAY_HANDLE_HPP
#define LNOT_ARRAY_HANDLE_HPP

#include <utility>

#include <LNOT/Containers/HostAllocator.hpp>
#include <LNOT/Settings.hpp>

namespace lnot
{
// NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic) : This class needs pointer arithmetic

/**
 * @brief Handle for managing dynamically allocated arrays using a specified allocator.
 * @tparam AllocatorT The type of the allocator to use for memory management.
 */
template <class AllocatorT>
class ArrayHandle
{
public:
  using Allocator = AllocatorT;                ///< The type of the allocator used for memory management.
  using Scalar = typename Allocator::Scalar;   ///< The type of the elements in the array.
  using Pointer = typename Allocator::Pointer; ///< Pointer to the allocated memory.
  using Size = typename Allocator::Size;       ///< The type used for sizes and indices.

  constexpr ArrayHandle() noexcept = default;
  constexpr ArrayHandle(const ArrayHandle&) noexcept = delete;
  constexpr ArrayHandle(ArrayHandle&& p_other) noexcept;

  /**
   * @brief Destructor for the ArrayHandle.
   * This will free the allocated memory if the handle is managing any.
   */
  ~ArrayHandle() noexcept { free(); }

  /**
   * @brief Delete copy assignment operator.
   * ArrayHandle cannot be copy-assigned to prevent multiple ownership of the same memory.
   */
  constexpr ArrayHandle& operator=(const ArrayHandle&) noexcept = delete;

  /**
   * @brief Move-assign another ArrayHandle to this one.
   * This will transfer ownership of the allocated memory from the other handle to this one.
   * The current handle will be freed if it is managing any memory.
   * After the move, the other handle will be set to a null state.
   * @param p_other The other ArrayHandle to move from.
   * @return Reference to this ArrayHandle.
   */
  constexpr ArrayHandle& operator=(ArrayHandle&& p_other) noexcept;

  /**
   * @brief Reset the ArrayHandle to a null state.
   * This will free the currently allocated memory and set the handle to nullptr.
   */
  constexpr ArrayHandle& operator=(std::nullptr_t /* p_handle */) noexcept;

  /**
   * @brief Checks if the ArrayHandle is managing a valid array.
   * @return True if the handle is not null, false otherwise.
   */
  constexpr explicit operator bool() const noexcept { return m_handle != nullptr; }

  /**
   * @brief Get the pointer to the allocated memory.
   * @return Pointer to the allocated memory.
   */
  constexpr Pointer cPtr() const noexcept { return m_handle; }

  /**
   * @brief Get the pointer to the allocated memory.
   * @return Pointer to the allocated memory.
   */
  constexpr Pointer operator->() const noexcept { return m_handle; }

  /**
   * @brief Dereference the pointer to the allocated memory.
   * @return Reference to the element pointed to by the handle.
   */
  constexpr Scalar& operator*() const noexcept { return *m_handle; }

  /**
   * @brief Access an element of the array by index.
   * @param i The index of the element to access.
   * @return Reference to the element at the specified index.
   */
  constexpr Scalar& operator[](Size i) const noexcept { return m_handle[i]; }

  /**
   * @brief Reallocate the array to hold n elements.
   * @param n The number of elements to allocate.
   */
  void reallocate(Size n) noexcept;

  /**
   * @brief Free the allocated memory.
   */
  void free() noexcept;

  /**
   * @brief Swap the handles of two ArrayHandle objects.
   * @param other The other ArrayHandle to swap with.
   */
  constexpr void swap(ArrayHandle& other) noexcept { std::swap(m_handle, other.m_handle); }

  /**
   * @brief Allocate a new ArrayHandle with the specified number of elements.
   * @param n The number of elements to allocate.
   * @return A new ArrayHandle managing the allocated memory.
   */
  static ArrayHandle allocate(Size n) noexcept { return ArrayHandle(Allocator::allocate(std::move(n))); }

private:
  explicit constexpr ArrayHandle(Pointer p_handle) noexcept;

  Pointer m_handle = nullptr;
};

// NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)

extern template class ArrayHandle<HostAllocator<settings::DefaultScalar, settings::DefaultSize>>;

/**
 * @brief Alias for an ArrayHandle using the HostAllocator with specified scalar and size types.
 */
template <typename ScalarT, typename SizeT = settings::DefaultSize>
using HostArrayHandle = ArrayHandle<HostAllocator<ScalarT, SizeT>>;

} // namespace lnot

#include <LNOT/Containers/ArrayHandle_impl.hpp>

#endif // LNOT_ARRAY_HANDLE_HPP
