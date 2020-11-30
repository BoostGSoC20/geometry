// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_FPG_ERROR_BOUND_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_FPG_ERROR_BOUND_HPP

#include <limits>
#include <cassert>
#include <array>
#include <tuple>
#include <algorithm>

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/set.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/map.hpp>
#include <boost/mp11/bind.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expansion_arithmetic.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/semi_static_filter.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template<typename ...Ts> using fpg_groups = std::tuple<Ts...>;
template<std::size_t ...> struct fpg_group {};

namespace fpg
{

// The following is an overestimation of ulp by at most a factor of 2.
// This could be improved.
template <typename Real>
constexpr Real ulp(Real d = Real(1))
{
    assert( d >= 0 );
    return d * std::numeric_limits<Real>::epsilon();
}

template <typename Real>
constexpr Real round_up_1_n(int n)
{
    Real out(1);
    for(int i = 0; i < n; ++i)
    {
        out *= (Real(1) + ulp<Real>());
    }
    return out;
}

template
<
    typename Real
>
struct static_filter_error
{
    Real magnitude;
    Real error;
};

template <typename Expression, bool = Expression::is_leaf>
constexpr bool is_arg = false;

template <typename Expression>
constexpr bool is_arg<Expression, true> = Expression::argn > 0;

template <typename Real>
constexpr Real sum_round_to_inf(Real a, Real b)
{
    assert( a >= 0 && b >= 0 );
    Real sum = a + b;
    Real tail = two_sum_tail(a, b, sum);
    if(tail > 0)
    {
        sum += ulp(sum);
    }
    return sum;
}

template <typename Real>
constexpr Real product_round_to_inf(Real a, Real b)
{
    assert( a >= 0 && b >= 0 );
    Real product = a * b;
    Real tail = two_product_tail_constexpr(a, b, product);
    if(tail > 0)
    {
        product += ulp(product);
    }
    return product;
}

template
<
    typename Expression,
    typename Real,
    operator_arities = Expression::operator_arity,
    bool IsArg = is_arg<Expression>
>
struct compute_static_filter_error;

template <typename Expression, typename Real>
struct compute_static_filter_error<Expression, Real, operator_arities::nullary, false>
{
    static constexpr static_filter_error<Real> apply()
    {
        return { std::abs(Expression::value), 0 };
    }
};

template <typename Expression, typename Real>
struct compute_static_filter_error<Expression, Real, operator_arities::nullary, true>
{
    static constexpr static_filter_error<Real> apply()
    {
        return { 1, 0 };
    }
};

template <typename Expression, typename Real>
struct compute_static_filter_error<Expression, Real, operator_arities::binary, false>
{
    static constexpr static_filter_error<Real> apply()
    {
        using l = typename Expression::left;
        using r = typename Expression::right;
        auto op = Expression::operator_type;
        if(op == operator_types::difference && is_arg<l> && is_arg<r>)
        {
            return { 1.0, ulp(Real(1)) / 2 };
        }
        else {
            constexpr auto e1 = compute_static_filter_error<l, Real>::apply();
            constexpr auto e2 = compute_static_filter_error<r, Real>::apply();
            if(op == operator_types::sum || op == operator_types::difference)
            {
                Real m = sum_round_to_inf(e1.magnitude, e2.magnitude);
                Real u = ulp(m) / 2;
                m = sum_round_to_inf(m, u);
                Real error = sum_round_to_inf(u, sum_round_to_inf(e1.error,
                                                                  e2.error));
                return { m, error };
            }
            else if (op == operator_types::product)
            {
                Real m = product_round_to_inf(e1.magnitude, e2.magnitude);
                Real u = ulp(m) / 2;
                m = sum_round_to_inf(m, u);
                Real error = sum_round_to_inf(u, product_round_to_inf(e1.error, e2.error));
                error = sum_round_to_inf(error, product_round_to_inf(e1.error, e2.magnitude));
                error = sum_round_to_inf(error, product_round_to_inf(e1.magnitude, e2.error));
                return { m, error };
            }
        }
        assert(false);
    }
};

constexpr int nonhomogenous = -1;

enum class decomposition_cases { general_binary, arg_diff, arg, constant, unhandled };

template
<
    typename Expression,
    operator_arities Arity = Expression::operator_arity
>
constexpr decomposition_cases decomposition_case = decomposition_cases::unhandled;

template <typename Expression>
constexpr decomposition_cases decomposition_case
    <
        Expression,
        operator_arities::binary
    > =    Expression::operator_type == operator_types::difference
        && is_arg<typename Expression::left>
        && is_arg<typename Expression::right> ?
          decomposition_cases::arg_diff
        : decomposition_cases::general_binary;

template <typename Expression>
constexpr decomposition_cases decomposition_case
    <
        Expression,
        operator_arities::nullary
    > = Expression::argn > 0 ? decomposition_cases::arg : decomposition_cases::constant;

template <typename Expression, decomposition_cases DC = decomposition_case<Expression>>
constexpr int degree = 1;

template <typename Expression>
constexpr int degree<Expression, decomposition_cases::constant> = 0;

template <typename Expression>
constexpr int degree<Expression, decomposition_cases::general_binary> =
    Expression::operator_type == operator_types::product ?
      degree<typename Expression::left> + degree<typename Expression::right>
    : (degree<typename Expression::left> == degree<typename Expression::right> ?
         degree<typename Expression::left>
       : nonhomogenous);

using arg_or_argdiff = std::array<std::size_t, 2>;

template <typename Expression, decomposition_cases = decomposition_case<Expression>>
constexpr int summands = 1;

template <typename Expression>
constexpr int summands<Expression, decomposition_cases::general_binary> =
    Expression::operator_type == operator_types::product ?
          summands<typename Expression::left> * summands<typename Expression::right>
        : summands<typename Expression::left> + summands<typename Expression::right>;

template
<
    typename Expression,
    std::size_t MaxArg,
    decomposition_cases = decomposition_case<Expression>
>
struct translation_groups_rec
{
    static constexpr std::array<std::size_t, MaxArg>
        apply(const std::array<std::size_t, MaxArg>& in)
    {
        return in;
    }
};

template <typename Expression, std::size_t MaxArg>
struct translation_groups_rec
    <
        Expression,
        MaxArg,
        decomposition_cases::general_binary
    >
{
    static constexpr std::array<std::size_t, MaxArg>
        apply(const std::array<std::size_t, MaxArg>& groups)
    {
        auto l_groups =
            translation_groups_rec<typename Expression::left, MaxArg>::apply(groups);
        return translation_groups_rec<typename Expression::right, MaxArg>
            ::apply(l_groups);
    }
};

template <std::size_t MaxArg, std::size_t ...Is>
static constexpr std::array<std::size_t, MaxArg>
    merge_groups(const std::array<std::size_t, MaxArg>& groups,
                 const std::size_t lower,
                 const std::size_t upper,
                 const std::index_sequence<Is...>)
{
    std::size_t out[MaxArg]{ (groups[Is])... };
    for(std::size_t i = 0; i < MaxArg; ++i)
    {
        if(out[i] == upper)
        {
            out[i] = lower;
        }
    }
    return std::array<std::size_t, MaxArg>{ (out[Is])... };
}

template <typename Expression, std::size_t MaxArg>
struct translation_groups_rec
    <
        Expression,
        MaxArg,
        decomposition_cases::arg_diff
    >
{
    static constexpr std::array<std::size_t, MaxArg>
        apply(const std::array<std::size_t, MaxArg>& groups)
    {
        std::size_t lower = std::min(groups[Expression::left::argn - 1],
                                     groups[Expression::right::argn - 1]);
        std::size_t upper = std::max(groups[Expression::left::argn - 1],
                                     groups[Expression::right::argn - 1]);
        if(upper > lower)
        {
            return merge_groups<MaxArg>(groups,
                                        lower,
                                        upper,
                                        std::make_index_sequence<MaxArg>{});
        }
        else
        {
            return groups;
        }
    }
};

template <std::size_t MaxArg, std::size_t ...Is>
constexpr std::array<std::size_t, MaxArg>
reindex_groups(const std::array<std::size_t, MaxArg>& groups,
               const std::index_sequence<Is...>)
{
    std::size_t new_index = 0;
    std::size_t cur = 1;
    std::size_t next = 0;
    std::size_t out[MaxArg] { (groups[Is])... };
    do
    {
        next = 0;
        for(std::size_t i = 0; i < MaxArg; ++i)
        {
            if(out[i] == cur)
            {
                out[i] = new_index;
            }
            else if(out[i] > cur && (next == 0 || next > out[i]))
            {
                next = out[i];
            }
        }
        cur = next;
        ++new_index;
    } while(next != 0);
    return std::array<std::size_t, MaxArg>{ (out[Is])... };
}

template <std::size_t MaxArg>
constexpr std::size_t
count_group_index_occurences(const std::array<std::size_t, MaxArg>& groups,
                             std::size_t group_index)
{
    std::size_t occurences = 0;
    for(std::size_t i = 0; i < MaxArg; ++i)
    {
        if(groups[i] == group_index)
        {
            ++occurences;
        }
    }
    return occurences;
}

template <std::size_t MaxArg, std::size_t ...Is>
constexpr std::array<std::size_t, MaxArg>
iota(const std::size_t offset, const std::index_sequence<Is...>)
{
    return std::array<std::size_t, MaxArg>{ (Is + offset)... };
}

template <std::size_t MaxArg>
constexpr std::size_t array_max(const std::array<std::size_t, MaxArg>& in)
{
    std::size_t m = 0;
    for(std::size_t i = 0; i < MaxArg; ++i)
    {
        m = std::max(m, in[i]);
    }
    return m;
}

template <typename Expression>
constexpr auto translation_group_assignments()
{
    constexpr std::size_t max_arg = max_argn<Expression>;
    constexpr std::make_index_sequence<max_arg> is{};
    constexpr auto groups_init = iota<max_arg>(1, is);
    constexpr auto groups = translation_groups_rec<Expression, max_arg>::apply(groups_init);
    constexpr auto groups_reindexed = reindex_groups<max_arg>(groups, is);
    return groups_reindexed;//groups_to_type(groups_reindexed, std::make_index_sequence<groups_count>{});
}

template <std::size_t GroupLength, std::size_t MaxArg, std::size_t ...Is>
constexpr std::array<std::size_t , GroupLength>
filter_group(const std::array<std::size_t, MaxArg> groups,
             const std::size_t group_index,
             std::index_sequence<Is...>)
{
    std::size_t out[GroupLength]{};
    std::size_t j = 0;
    for(std::size_t i = 0; i < MaxArg; ++i)
    {
        if(groups[i] == group_index)
        {
            out[j++] = i + 1;
        }
    }
    return std::array<std::size_t, GroupLength>{ (out[Is])... };
}

template <std::size_t Length, std::size_t ...Is>
constexpr auto array_to_fpg_group(const std::array<std::size_t, Length> in,
                                  const std::index_sequence<Is...>)
{
    return fpg_group< (in[Is])... >{};
}

template <typename Expression, std::size_t GroupIndex>
constexpr auto translation_group_args()
{
    constexpr auto groups = translation_group_assignments<Expression>();
    constexpr std::size_t group_length =
        count_group_index_occurences(groups, GroupIndex);
    using group_indices = std::make_index_sequence<group_length>;
    return filter_group<group_length>(groups, GroupIndex, group_indices{});
}

template <typename Expression, std::size_t GroupIndex, std::size_t ...Is>
constexpr auto translation_fpg_group(std::index_sequence<Is...>)
{
    constexpr auto args = translation_group_args<Expression, GroupIndex>();
    return fpg_group< (args[Is])... >{};
}

template <typename Expression, std::size_t GroupIndex>
constexpr std::size_t translation_group_length =
    count_group_index_occurences(translation_group_assignments<Expression>(), GroupIndex);

template <typename Expression, std::size_t GroupIndex>
using translation_fpg_group_t =
    decltype(translation_fpg_group<Expression, GroupIndex>(
        std::make_index_sequence<translation_group_length<Expression, GroupIndex>>{}));

template <typename Expression>
constexpr std::size_t translation_groups_count =
    array_max(translation_group_assignments<Expression>()) + 1;

template
<
    typename Expression,
    std::size_t I = translation_groups_count<Expression>,
    typename ...FpgGroups
>
struct translation_fpg_groups_impl
{
    using type = typename translation_fpg_groups_impl
        <
            Expression,
            I - 1,
            translation_fpg_group_t<Expression, I - 1>,
            FpgGroups...
        >::type;
};

template
<
    typename Expression,
    typename ...FpgGroups
>
struct translation_fpg_groups_impl<Expression, 0, FpgGroups...>
{
    using type = fpg_groups<FpgGroups...>;
};

template <typename Expression>
using translation_fpg_groups_t = typename translation_fpg_groups_impl<Expression>::type;

template <typename T, std::size_t L1, std::size_t L2, std::size_t ...Is>
constexpr std::array<T, L1 + L2> concat_arrays(const std::array<T, L1> arr1,
                                               const std::array<T, L2> arr2,
                                               const std::index_sequence<Is...>)
{
    T out[L1 + L2]{};
    std::size_t i_out = 0;
    for(std::size_t i = 0; i < L1; ++i)
    {
        out[i_out++] = arr1[i];
    }
    for(std::size_t i = 0; i < L2; ++i)
    {
        out[i_out++] = arr2[i];
    }
    return std::array<T, L1 + L2>{ (out[Is])... };
}

template <typename T, std::size_t LI1, std::size_t LI2, std::size_t LO1, std::size_t LO2, std::size_t ...Is>
constexpr std::array<std::array<T, LI1 + LI2>, LO1 * LO2>
array_array_product(const std::array<std::array<T, LI1>, LO1> arr1,
                    const std::array<std::array<T, LI2>, LO2> arr2,
                    const std::index_sequence<Is...>)
{
    return std::array<std::array<T, LI1 + LI2>, LO1 * LO2>
        { (concat_arrays(arr1[Is % LO1], arr2[Is / LO1], std::make_index_sequence<LI1 + LI2>{}))... };
}

template
<
    typename Expression,
    decomposition_cases = decomposition_case<Expression>,
    operator_types = Expression::operator_type
>
struct expand_polynomial_impl;

template
<
    typename Expression,
    operator_types Op
>
struct expand_polynomial_impl<Expression, decomposition_cases::arg, Op>
{
    static constexpr std::array<std::array<arg_or_argdiff, 1>, 1> apply()
    {
        return std::array<std::array<arg_or_argdiff, 1>, 1>{
            std::array<arg_or_argdiff, 1>{ arg_or_argdiff{Expression::argn, 0} }
        };
    }
};

template <typename Expression, operator_types Op>
struct expand_polynomial_impl<Expression, decomposition_cases::constant, Op>
{
    static constexpr std::array<std::array<arg_or_argdiff, 0>, 1> apply()
    {
        return std::array<std::array<arg_or_argdiff, 0>, 1>{ std::array<arg_or_argdiff, 0>{} };
    }
};

template <typename Expression, operator_types Op>
struct expand_polynomial_impl<Expression, decomposition_cases::arg_diff, Op>
{
    static constexpr std::array<std::array<arg_or_argdiff, 1>, 1> apply()
    {
        return std::array<std::array<arg_or_argdiff, 1>, 1>{
            std::array<arg_or_argdiff, 1>{
                arg_or_argdiff{
                    Expression::left::argn,
                    Expression::right::argn
                }
            }
        };
    }
};

template <typename Expression>
struct expand_polynomial_impl
    <
        Expression,
        decomposition_cases::general_binary,
        operator_types::product
    >
{
    static constexpr std::array<
            std::array<arg_or_argdiff, degree<Expression>>,
            summands<Expression>
        > apply()
        {
            constexpr auto l_expanded = expand_polynomial_impl<typename Expression::left>::apply();
            constexpr auto r_expanded = expand_polynomial_impl<typename Expression::right>::apply();
            return array_array_product(
                    l_expanded,
                    r_expanded,
                    std::make_index_sequence<l_expanded.size() * r_expanded.size()>{});
        }
};

template <typename Expression, operator_types Op>
struct expand_polynomial_impl<Expression, decomposition_cases::general_binary, Op>
{
    static constexpr std::array<
            std::array<arg_or_argdiff, degree<Expression>>,
            summands<Expression>
        > apply()
        {
            constexpr auto l_expanded = expand_polynomial_impl<typename Expression::left>::apply();
            constexpr auto r_expanded = expand_polynomial_impl<typename Expression::right>::apply();
            return concat_arrays(
                    l_expanded,
                    r_expanded,
                    std::make_index_sequence<l_expanded.size() + r_expanded.size()>{});
        }
};

template <std::size_t ...GroupArgs>
constexpr bool is_in_group(const arg_or_argdiff& argd, const fpg_group<GroupArgs...>)
{
    bool argd1 = false;
    bool argd2 = argd[1] != 0;
    for(const auto group_arg : { GroupArgs... })
    {
        argd1 = argd1 || argd[0] == group_arg;
        argd2 = argd2 || argd[1] == group_arg;
    }
    return argd1 && argd2;
}

template
<
    std::size_t Summands,
    std::size_t Degree,
    std::size_t ...GroupArgs
>
constexpr void assign_to_group(
    const std::array<std::array<arg_or_argdiff, Degree>, Summands> expanded_polynomial,
    std::size_t group_assignments[Summands][Degree],
    const fpg_group<GroupArgs...> group,
    const std::size_t group_index)
{
    bool group_factor_remaining = true;
    do
    {
        group_factor_remaining = true;
        bool summand_without_group_factor = false;
        for(std::size_t i = 0; i < Summands; ++i)
        {
            bool group_factor_in_summand = false;
            for(std::size_t j = 0; j < Degree; ++j)
            {
                if(   group_assignments[i][j] == 0
                   && is_in_group(expanded_polynomial[i][j], group))
                {
                    group_factor_in_summand = true;
                    break;
                }
            }
            if(!group_factor_in_summand)
            {
                group_factor_remaining = false;
                break;
            }
        }
        if( !group_factor_remaining )
        {
            break;
        }
        for(std::size_t i = 0; i < Summands; ++i)
        {
            for(std::size_t j = 0; j < Degree; ++j)
            {
                if(   group_assignments[i][j] == 0
                   && is_in_group(expanded_polynomial[i][j], group))
                {
                    group_assignments[i][j] = group_index + 1;
                    break;
                }
            }
        }
    } while(group_factor_remaining);
}

template
<
    std::size_t Summands,
    std::size_t Degree,
    typename ...FpgGroups,
    std::size_t ...Is, template <typename, std::size_t ...> class T,
    std::size_t ...Js, template <typename, std::size_t ...> class V
>
constexpr std::array<std::array<std::size_t, Degree>, Summands> assign_to_groups(
    const std::array<std::array<arg_or_argdiff, Degree>, Summands>& expanded_polynomial,
    const fpg_groups<FpgGroups...>,
    const T<std::size_t, Is...>,
    const V<std::size_t, Js...>)
{
    static_assert(sizeof...(Is) == sizeof...(FpgGroups), "Mismatching index sequence");
    static_assert(sizeof...(Js) == Summands * Degree, "Mismatching index sequence");
    std::size_t group_assignments[Summands][Degree]{};
    for(std::size_t i = 0; i < Summands; ++i)
    {
        for(std::size_t j = 0; j < Degree; ++j)
        {
            group_assignments[i][j] = 0;
        }
    }
    using dummy = int[];
    (void)dummy {0, ((void)assign_to_group(expanded_polynomial, group_assignments, FpgGroups{}, Is), 0)...};
    return std::array<std::array<std::size_t, Degree>, Summands>
        { (group_assignments[Js / Degree][Js % Degree])... };
}

template <std::size_t Summands, std::size_t Degree>
constexpr std::size_t group_degree(
    const std::array<std::array<std::size_t, Degree>, Summands>& group_assignments,
    const std::size_t group_index)
{
    std::size_t deg = 0;
    for(std::size_t i = 0; i < Degree; ++i)
    {
        if(group_assignments[0][i] == group_index)
        {
            ++deg;
        }
    }
    return deg;
}

template <std::size_t OutL, std::size_t Degree, std::size_t Summands, std::size_t ...Is>
constexpr std::array<arg_or_argdiff, OutL>
extract_by_group(const std::array<std::array<arg_or_argdiff, Degree>, Summands>& expanded_polynomial,
                 const std::array<std::array<std::size_t, Degree>, Summands>& group_assignments,
                 const std::size_t group_index,
                 const std::index_sequence<Is...>)
{
    static_assert(OutL == sizeof...(Is), "Index sequence must match output length");
    std::size_t out[OutL][2]{};
    std::size_t i_out = 0;
    for(int i = 0; i < Summands; ++i)
    {
        for(int j = 0; j < Degree; ++j)
        {
            if(group_assignments[i][j] == group_index)
            {
                out[i_out][0] = expanded_polynomial[i][j][0];
                out[i_out++][1] = expanded_polynomial[i][j][1];
            }
        }
    }
    return std::array<arg_or_argdiff, OutL>
        { (arg_or_argdiff{out[Is][0], out[Is][1]})... };
}

template <std::size_t L>
constexpr std::size_t count_unique(const std::array<arg_or_argdiff, L>& arr)
{
    std::size_t uniques = 0;
    for(std::size_t i = 0; i < arr.size(); ++i)
    {
        bool unique = true;
        for(std::size_t j = 0; j < i; ++j)
        {
            if(arr[i][0] == arr[j][0] && arr[i][1] == arr[j][1])
            {
                unique = false;
            }
        }
        if(unique)
        {
            ++uniques;
        }
    }
    return uniques;
}

template
<
    std::size_t U,
    std::size_t L,
    std::size_t ...Is
>
constexpr std::array<arg_or_argdiff, U>
filter_unique(const std::array<arg_or_argdiff, L>& arr,
              const std::index_sequence<Is...>)
{
    static_assert(U == sizeof...(Is), "Index sequence must match output length.");
    std::size_t out[U][2]{};
    std::size_t i_out = 0;
    for(std::size_t i = 0; i < arr.size(); ++i)
    {
        bool unique = true;
        for(std::size_t j = 0; j < i; ++j)
        {
            if(arr[i][0] == arr[j][0] && arr[i][1] == arr[j][1])
            {
                unique = false;
            }
        }
        if(unique)
        {
            out[i_out][0] = arr[i][0];
            out[i_out++][1] = arr[i][1];
        }
    }
    return std::array<arg_or_argdiff, U>
        { (arg_or_argdiff{out[Is][0], out[Is][1]})... };
}

template <typename Expression, typename FpgGroups, std::size_t GroupIndex>
constexpr auto group_arg_or_argdiffs()
{
    constexpr auto expanded = fpg::expand_polynomial_impl<Expression>::apply();
    constexpr auto group_assignments =
        fpg::assign_to_groups(
            expanded,
            FpgGroups{},
            std::make_index_sequence<std::tuple_size<FpgGroups>::value>{},
            std::make_index_sequence<expanded.size() * expanded[0].size()>{});
    constexpr std::size_t gd = group_degree(group_assignments, GroupIndex);
    constexpr auto all_group_factors =
        extract_by_group<gd * expanded.size()>(expanded,
                                               group_assignments,
                                               GroupIndex,
                                               std::make_index_sequence<gd * expanded.size()>{});
    constexpr std::size_t unique_factors_count = count_unique(all_group_factors);
    constexpr auto unique_factors =
        filter_unique<unique_factors_count>(all_group_factors,
                                            std::make_index_sequence<unique_factors_count>{});
    return std::make_pair(gd, unique_factors);
}

template <std::size_t argn1, std::size_t argn2>
struct arg_or_arg_diff_exp_impl
{
    using type = abs<difference< argument<argn1>, argument<argn2> >>;
};

template <std::size_t argn>
struct arg_or_arg_diff_exp_impl<argn, 0>
{
    using type = abs<argument<argn>>;
};

template <typename ...Exps>
struct multi_max_impl {};

template <typename Exp1, typename ...Exps>
struct multi_max_impl<Exp1, Exps...>
{
    using type = max<Exp1, typename multi_max_impl<Exps...>::type>;
};

template <typename Exp>
struct multi_max_impl<Exp>
{
    using type = Exp;
};

template <>
struct multi_max_impl<>
{
    using type = int;
};

template <typename Exp, std::size_t exp>
struct power_impl
{
    using type = product<Exp, typename power_impl<Exp, exp - 1>::type>;
};

template <typename Exp>
struct power_impl<Exp, 1>
{
    using type = Exp;
};

template <typename Exp>
struct power_impl<Exp, 0>
{
    using type = int;
};

template <typename Expression, typename FpgGroups, std::size_t GroupIndex, std::size_t ...Is>
constexpr auto group_bound_expression_helper(const std::index_sequence<Is...>)
{
    constexpr auto group_info =
        group_arg_or_argdiffs<Expression, FpgGroups, GroupIndex>();
    static_assert(sizeof...(Is) == group_info.second.size(), "Index sequence length must match group factor count.");
    using mm = typename multi_max_impl
        <
            typename arg_or_arg_diff_exp_impl
                <
                    group_info.second[Is][0],
                    group_info.second[Is][1]
                >::type...
        >::type;
    return typename power_impl<mm, group_info.first>::type{};
}

template <typename Expression, typename FpgGroups, std::size_t GroupIndex>
constexpr auto group_bound_expression()
{
    constexpr auto group_info =
        group_arg_or_argdiffs<Expression, FpgGroups, GroupIndex>();
    return group_bound_expression_helper<Expression, FpgGroups, GroupIndex>(
        std::make_index_sequence<group_info.second.size()>{});
}

template <typename Expression, typename FpgGroups, std::size_t GroupIndex>
using group_bound_expression_t =
    decltype(group_bound_expression<Expression, FpgGroups, GroupIndex>());

template <typename Expression, typename FpgGroups, std::size_t ...Is>
constexpr auto groups_bound_expression_helper(const std::index_sequence<Is...>)
{
    using all_bounds = boost::mp11::mp_remove
        <
            boost::mp11::mp_list<group_bound_expression_t<Expression, FpgGroups, Is>...>,
            int
        >;
    using first_bound = boost::mp11::mp_front<all_bounds>;
    using remainder = boost::mp11::mp_pop_front<all_bounds>;
    using result = boost::mp11::mp_fold<remainder, first_bound, product>;
    return result{};
}

template <typename Expression, typename FpgGroups>
using groups_bound_expression_t =
    decltype(groups_bound_expression_helper<Expression, FpgGroups>(
        std::make_index_sequence<boost::mp11::mp_size<FpgGroups>::value + 1>{}));


// The following template derives an error expression inspired by the ideas of
// "FPG: A code generator for fast and certified geometric predicates" by Meyer
// and Pion. The implementation (at the time of writing this comment) makes the
// following assumptions:
// 1. Groups is an exact cover of the set of arguments contained in Expression
// 2. The expanded polynomial consists of summands such that for each group
//    each summand contains the same number of factors out of that group.
// 3. There are no higher-degree arguments.
template
<
    typename Expression,
    typename Real,
    typename Groups = translation_fpg_groups_t<Expression>
>
struct error_expression
{
private:
    static constexpr std::size_t deg = degree<Expression>;
    static constexpr Real delta_1 =
          compute_static_filter_error<Expression, Real>::apply().error
        * round_up_1_n<Real>(deg);
    struct delta_constant : public static_constant_interface<Real>
    {
        static constexpr Real value = delta_1;
        static constexpr bool non_negative = true;
    };
    using scale = groups_bound_expression_t<Expression, Groups>;
public:
    using type = product<delta_constant, scale>;
};

} // fpg

template
<
    typename Expression,
    typename CalculationType,
    typename Groups = fpg::translation_fpg_groups_t<Expression>
>
using fpg_error_expression = typename fpg::error_expression
        <
            Expression,
            CalculationType,
            Groups
        >::type;

template
<
    typename Expression,
    typename CalculationType,
    typename Groups = fpg::translation_fpg_groups_t<Expression>
>
using fpg_semi_static = semi_static_filter
        <
            Expression,
            CalculationType,
            fpg_error_expression
                <
                    Expression,
                    CalculationType,
                    Groups
                >
        >;
}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_FPG_ERROR_BOUND_HPP
