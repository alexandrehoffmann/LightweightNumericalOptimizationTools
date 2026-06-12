#ifndef LNOT_ASSERT_HPP
#define LNOT_ASSERT_HPP

#include <cstdio>
#include <cstdlib>

namespace LNOT
{

inline void lnot_assert(const bool cond, const char* msg, const char* file, const int line)
{
	if (!cond) 
	{
		std::fprintf(stderr, "LNOT assertion failed: %s\n  %s:%d\n", msg, file, line);
		std::abort();
	}
}

#define LNOT_ASSERT(cond, msg) lnot_assert(cond, msg, __FILE__, __LINE__)

} // namespace LNOT

#endif // LNOT_ASSERT_HPP
