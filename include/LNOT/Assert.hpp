#ifndef LNOT_ASSERT_HPP
#define LNOT_ASSERT_HPP

#include <exception>
#include <print>
#include <string_view>

namespace lnot
{

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters,-warnings-as-errors)
constexpr void lnotAssert(const bool p_cond, std::string_view p_msg, std::string_view p_file = __FILE__, const int p_line = __LINE__)
{
  if (!p_cond)
  {
    std::println(stderr, "LNOT assertion failed: {}\n  {}:{}", p_msg, p_file, p_line);
    std::terminate();
  }
}

} // namespace lnot

#endif // LNOT_ASSERT_HPP
