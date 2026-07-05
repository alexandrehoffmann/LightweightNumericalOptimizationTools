#ifndef LNOT_TRAITS_HPP
#define LNOT_TRAITS_HPP

#include <type_traits>

namespace LNOT
{

namespace detail
{

template<typename T> concept CWithScalar    = requires { typename std::decay_t<T>::Scalar;    };
template<typename T> concept CWithSize      = requires { typename std::decay_t<T>::Size;      };
template<typename T> concept CWithCriterion = requires { typename std::decay_t<T>::Criterion; };
template<typename T> concept CWithInfo      = requires { typename std::decay_t<T>::Info;      };

} // namespace detail

template<detail::CWithScalar T>    using ScalarFor    = typename std::decay_t<T>::Scalar;
template<detail::CWithSize T>      using SizeFor      = typename std::decay_t<T>::Size;
template<detail::CWithCriterion T> using CriterionFor = typename std::decay_t<T>::Criterion;
template<detail::CWithCriterion T> using InfoFor      = typename std::decay_t<T>::Info;

} // namespace LNOT

#endif // LNOT_TRAITS_HPP
