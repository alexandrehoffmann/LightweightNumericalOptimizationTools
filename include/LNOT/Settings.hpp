#ifndef LNOT_SETTINGS_HPP
#define LNOT_SETTINGS_HPP

namespace LNOT::settings
{

#ifdef LNOT_WITH_BLAS
inline constexpr bool useBlas = true;
#else
inline constexpr bool useBlas = false;
#endif // LNOT_WITH_BLAS

#ifdef LNOT_WITH_CUDA
inline constexpr bool useCuda = true;
#else
inline constexpr bool useCuda = false;
#endif // LNOT_WITH_BLAS

} // namespace LNOT::settings

#endif // LNOT_SETTINGS_HPP
