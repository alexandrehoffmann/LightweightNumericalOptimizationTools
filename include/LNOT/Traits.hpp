#ifndef LNOT_TRAITS_HPP
#define LNOT_TRAITS_HPP

namespace LNOT
{

namespace detail
{

template<typename T> concept CWithScalar    = requires { typename T::Scalar;    };
template<typename T> concept CWithSize      = requires { typename T::Size;      };
template<typename T> concept CWithCriterion = requires { typename T::Criterion; };
template<typename T> concept CWithInfo      = requires { typename T::Info;      };

} // namespace detail

template<detail::CWithScalar T>    using ScalarFor    = typename T::Scalar;
template<detail::CWithSize T>      using SizeFor      = typename T::Size;
template<detail::CWithCriterion T> using CriterionFor = typename T::Criterion;
template<detail::CWithCriterion T> using InfoFor      = typename T::Info;

} // namespace LNOT

#endif // LNOT_TRAITS_HPP
