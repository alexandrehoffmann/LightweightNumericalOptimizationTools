#ifndef LNOT_HOST_ALLOCATOR_HPP
#define LNOT_HOST_ALLOCATOR_HPP

#include <new>
#include <utility>

#include <LNOT/Settings.hpp>

namespace lnot
{

// NOLINTBEGIN(cppcoreguidelines-owning-memory) : Entire point of the class
/**
 * @brief Host memory allocator for LNOT containers.
 * @tparam ScalarT The type of the elements to allocate.
 * @tparam SizeT The type used for sizes and indices.
 */
template <typename ScalarT, typename SizeT = settings::DefaultSize>
struct HostAllocator
{
  using Scalar = ScalarT;  ///< The type of the elements to allocate.
  using Size = SizeT;      ///< The type used for sizes and indices.
  using Pointer = Scalar*; ///< Pointer to the allocated memory.

  /**
   * @brief Allocate memory for n elements.
   * @param n The number of elements to allocate.
   * @return Pointer to the allocated memory.
   */
  static Pointer allocate(Size p_size) noexcept { return new (std::nothrow) Scalar[std::move(p_size)]; }

  /**
   * @brief Free the allocated memory.
   * @param ptr Pointer to the memory to free.
   */
  static void free(Pointer p_ptr) noexcept { delete[] p_ptr; }
};
// NOLINTEND(cppcoreguidelines-owning-memory)

extern template struct HostAllocator<settings::DefaultScalar, settings::DefaultSize>;

} // namespace lnot

#endif // LNOT_HOST_ALLOCATOR_HPP
