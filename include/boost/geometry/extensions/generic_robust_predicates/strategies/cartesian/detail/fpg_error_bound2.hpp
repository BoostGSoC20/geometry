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

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/set.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/map.hpp>
#include <boost/mp11/bind.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/semi_static_filter.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template<typename ...> struct fpg_groups {};
template<typename ...> struct fpg_group {};

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
    Real tail = two_product_tail(a, b, product);
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
    operator_arities = Expression::operator_arity
>
struct compute_static_filter_error;

template <typename Expression, typename Real>
struct compute_static_filter_error<Expression, Real, operator_arities::nullary>
{
    static constexpr static_filter_error<Real> apply()
    {
        if(is_arg<Expression>)
        {
            return { 1.0, 0 };
        }
        else
        {
            return { std::abs(Expression::value), 0 };
        }
    }
};

template <typename Expression, typename Real>
struct compute_static_filter_error<Expression, Real, operator_arities::binary>
{
    static constexpr static_filter_error<Real> apply()
    {
        using l = typename Expression::left;
        using r = typename Expression::right;
        using op = typename Expression::operator_type;
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

template
<
    typename Expression,
    bool Translation,
    operator_types Op = Expression::operator_type
>
constexpr bool expansion_anchor = Expression::is_leaf;

template <typename Expression, true, operator_types::difference>
constexpr bool expansion_anchor =    is_arg<typename Expression::left>
                                  && is_arg<typename Expression::right>;

template
<
    typename Expression,
    bool Translation = true,
    operator_types Op = Expression::operator_type,
    typename Anchor = expansion_anchor<Expression, Translation>
>
struct expand_polynomial_impl {};

template <typename ...> struct expand_sum {};
template <typename ...> struct expand_product {};

template
<
    typename Expression,
    bool Translation,
    operator_types Op
>
struct expand_polynomial_impl
    <
        Expression,
        Translation,
        Op,
        boost::mp11::mp_true
    >
{
    using type = expand_sum<expand_product<Expression>>;
};

template
<
    typename Expression,
    bool Translation
>
struct expand_polynomial_impl
    <
        Expression,
        Translation,
        operator_types::product,
        boost::mp11::mp_false
    >
{
private:
    using ld = typename expand_polynomial_impl
        <
            typename Expression::left,
            Translation
        >::type;
    using rd = typename expand_polynomial_impl
        <
            typename Expression::right,
            Translation
        >::type;
public:
    using type = boost::mp11::mp_product<boost::mp11::mp_append, ld, rd>;
};

template
<
    typename Expression,                                   
    bool Translation
>
struct expand_polynomial_impl
    <
        Expression,
        Translation,
        operator_types::sum,
        boost::mp11::mp_false
    >
{
private:
    using ld = typename expand_polynomial_impl
        <
            typename Expression::left,
            Translation
        >::type;
    using rd = typename expand_polynomial_impl
        <
            typename Expression::right,
            Translation
        >::type;
public:
    using type = boost::mp11::mp_append<ld, rd>;
};

template
<
    typename Expression,
    bool Translation
>
struct expand_polynomial_impl
    <
        Expression,
        Translation,
        operator_types::difference,
        boost::mp11::mp_false
    >
{
private:
    using ld = typename expand_polynomial_impl
        <
            typename Expression::left,
            Translation
        >::type;
    using rd = typename expand_polynomial_impl
        <
            typename Expression::right,
            Translation
        >::type;
public:
    using type = boost::mp11::mp_append<ld, rd>;
};

template <typename Expression, bool Translation = true>
using expand_polynomial =
    typename expand_polynomial_impl<Expression, Translation>::type;

template
<
    typename Arg,
    typename DegreeMap,
    typename Contains = boost::mp11::mp_map_contains<DegreeMap, Arg>
>
struct degree_helper
{
    using type =
        boost::mp11::mp_second<boost::mp11::mp_map_find<DegreeMap, Arg>>;
};

template <typename Arg, typename DegreeMap>
struct degree_helper<Arg, DegreeMap, boost::mp11::mp_false>
{
    using type = boost::mp11::mp_int<1>;
};

template
<
    typename LeafOrDiff,
    typename DegreeMap,
    bool IsLeaf = LeafOrDiff::is_leaf
>
struct degree_impl
{
    using type = typename degree_helper<LeafOrDiff, DegreeMap>::type;
};

template <typename Diff, typename DegreeMap>
struct degree_impl<Diff, DegreeMap, false>
{
private:
    using left_degree = degree_helper<typename Diff::left, DegreeMap>;
    using right_degree = degree_helper<typename Diff::right, DegreeMap>;
    static_assert(left_degree::type::value == right_degree::type::value,
                  "Must not mix degrees in difference.");
public:
    using type = typename left_degree::type;
};

template <typename LeafOrDiff, typename DegreeMap>
using degree = typename degree_impl<LeafOrDiff, DegreeMap>::type;

template <typename ExpandProduct, typename DegreeMap>
struct product_degree_impl
{
private:
    using degree_q = boost::mp11::mp_bind_back<degree, DegreeMap>;
    using degrees = boost::mp11::mp_transform_q<degree_q, ExpandProduct>;
public:
    using type = boost::mp11::mp_fold
        <
            degrees,
            boost::mp11::mp_int<0>,
            boost::mp11::mp_plus
        >;
};

template <typename ExpandProduct, typename DegreeMap>
using product_degree =
    typename product_degree_impl<ExpandProduct, DegreeMap>::type;

template
<
    typename ExpandPolynomial,
    typename DegreeMap = boost::mp11::mp_list<>
>
struct expand_polynomial_degree_impl
{
public:
    using type = typename product_degree_impl
        <
            boost::mp11::mp_first<ExpandPolynomial>,
            DegreeMap
        >::type;
private:
    using degree_q = boost::mp11::mp_bind_back<product_degree, DegreeMap>;
    using degrees = boost::mp11::mp_transform_q<degree_q, ExpandPolynomial>;
    template <typename Degree>
    using correct = boost::mp11::mp_bool<Degree::value == type::value>;
    using all_correct = boost::mp11::mp_all_of<degrees, correct>;
    static_assert(all_correct::value, "Polynomial must be homogenous.");
};

template<typename ExpandPoly, typename DegreeMap = boost::mp11::mp_list<>>
using expand_polynomial_degree =
    typename expand_polynomial_degree_impl<ExpandPoly, DegreeMap>::type;

template <typename Expression>
using is_difference = boost::mp11::mp_bool
    <
        Expression::operator_type == operator_types::difference
    >;

template <typename ArgGroupIndexMap, std::size_t NextIndex = 1>
struct translation_group_fold_state
{
    using map = ArgGroupIndexMap;
    static constexpr std::size_t next = NextIndex;
};

template
<
    typename State,
    typename Diff,
    bool known_group = boost::mp11::mp_or
        <
            boost::mp11::mp_map_contains
                <
                    typename State::map,
                    boost::mp11::mp_front<Diff>
                >,
            boost::mp11::mp_map_contains
                <
                    typename State::map,
                    boost::mp11::mp_back<Diff>
                >
        >::value
>
struct translation_group_fold_operation_impl
{
private:
    using key = boost::mp11::mp_if
        <
            boost::mp11::mp_map_contains
                <
                    typename State::map,
                    boost::mp11::mp_front<Diff>
                >,
            boost::mp11::mp_front<Diff>,
            boost::mp11::mp_back<Diff>
        >;
    using value = boost::mp11::mp_second
        <
            boost::mp11::mp_map_find<typename State::map, key>
        >;
    using inserted1 = boost::mp11::mp_map_insert
        <
            typename State::map,
            boost::mp11::mp_list
                <
                    boost::mp11::mp_second<Diff>,
                    value
                >
        >;
    using inserted2 = boost::mp11::mp_map_insert
        <
            inserted1,
            boost::mp11::mp_list
                <
                    boost::mp11::mp_first<Diff>,
                    value
                >
        >;
public:
    using type = translation_group_fold_state
        <
            inserted2,
            State::next
        >;
};

template <typename State, typename Diff>
struct translation_group_fold_operation_impl<State, Diff, false>
{
private:
    using key1 = boost::mp11::mp_front<Diff>;
    using key2 = boost::mp11::mp_back<Diff>;
    using value = boost::mp11::mp_size_t<State::next>;
    using inserted1 = boost::mp11::mp_map_insert
        <
            typename State::map,
            boost::mp11::mp_list
                <
                    key1,
                    value
                >
        >;
    using inserted2 = boost::mp11::mp_map_insert
        <
            inserted1,
            boost::mp11::mp_list
                <
                    key2,
                    value
                >
        >;
public:
    using type = translation_group_fold_state
        <
            inserted2,
            State::next + 1
        >;
};

template <typename State, typename Diff>
using translation_group_fold_operation =
    typename translation_group_fold_operation_impl<State, Diff>::type;

template <typename Expression>
struct translation_auto_groups
{
private:
    using expanded_polynomial = expand_polynomial<Expression>;
    using translations =
        boost::mp11::mp_unique
            <
                boost::mp11::mp_copy_if
                    <
                        boost::mp11::mp_flatten
                            <
                                expanded_polynomial,
                                expand_product<>
                            >,
                        is_difference
                    >
            >;
    using arg_to_group_index_map = boost::mp11::mp_fold
        <
            translations,
            translation_group_fold_state
                <
                    boost::mp11::mp_list<>,
                    0
                >,
            translation_group_fold_operation
        >;
    static constexpr std::size_t group_count = arg_to_group_index_map::next;
    using empty_groups = boost::mp11::mp_repeat_c
        <
            fpg_groups<fpg_group<>>,
            group_count
        >;
    template<typename Groups, typename ArgIndexPair>
    using build_groups_fold = boost::mp11::mp_replace_at
        <
            Groups,
            boost::mp11::mp_second<ArgIndexPair>,
            boost::mp11::mp_push_back
                <
                    boost::mp11::mp_at
                        <
                            Groups,
                            boost::mp11::mp_second<ArgIndexPair>
                        >,
                    boost::mp11::mp_first<ArgIndexPair>
                >
        >;
public:
    using type = boost::mp11::mp_fold
        <
            typename arg_to_group_index_map::map,
            empty_groups,
            build_groups_fold
        >;
};

template <typename Expression, std::size_t Exponent>
struct power_c_impl
{
private:
    using remainder = boost::mp11::mp_repeat_c
        <
            boost::mp11::mp_list<Expression>,
            Exponent - 1
        >;
public:
    using type = boost::mp11::mp_fold
        <
            remainder,
            Expression,
            product
        >;
};

template <typename Expression, std::size_t Exponent>
using power_c = typename power_c_impl<Expression, Exponent>::type;

template <typename Expression, typename Exponent>
using power = power_c<Expression, Exponent::value>;

template
<
    typename ArgOrDiff,
    typename Group,
    bool IsDiff = ArgOrDiff::operator_type == operator_types::difference
>
struct is_in_group_impl
{
    using type = boost::mp11::mp_contains<Group, ArgOrDiff>;
};

template
<
    typename Diff,
    typename Group
>
struct is_in_group_impl<Diff, Group, true>
{
private:
    using first = boost::mp11::mp_first<Diff>;
    using second = boost::mp11::mp_second<Diff>;
public:
    using type = boost::mp11::mp_or
        <
            boost::mp11::mp_contains<Group, Diff>,
            boost::mp11::mp_and
                <
                    boost::mp11::mp_contains<Group, first>,
                    boost::mp11::mp_contains<Group, second>
                >
        >;
};

template <typename ArgOrDiff, typename Group>
using is_in_group = typename is_in_group_impl<ArgOrDiff, Group>::type;

template <typename Group>
struct is_in_group_helper
{
    template <typename ArgOrDiff>
    using fn = is_in_group<ArgOrDiff, Group>;
};

template
<
    typename ExpandProduct
>
struct group_degree_helper
{
    template <typename Group>
    using fn = boost::mp11::mp_count_if_q
        <
            ExpandProduct,
            is_in_group_helper<Group>
        >;
};

template <typename Group, typename ExpandPolynomial>
struct group_expression_impl
{
private:
    using all_leaves = boost::mp11::mp_unique
        <
            boost::mp11::mp_flatten<ExpandPolynomial, expand_product<>>
        >;
    using in_group = is_in_group_helper<Group>;
    using leaves = boost::mp11::mp_copy_if_q<all_leaves, in_group>;
    using aleaves = boost::mp11::mp_transform<abs, leaves>;
public:
    using type = boost::mp11::mp_fold
        < 
            boost::mp11::mp_pop_front<aleaves>,
            boost::mp11::mp_front<aleaves>,
            max
        >;
};

template <typename Group, typename ExpandPolynomial>
using group_expression =
    typename group_expression_impl<Group, ExpandPolynomial>::type;

template
<
    typename ExpandPolynomial
>
struct group_expression_helper
{
    template <typename Group>
    using fn = group_expression<Group, ExpandPolynomial>;
};

template <typename Groups>
struct group_degrees_helper
{
    template <typename Summand>
    using fn = boost::mp11::mp_transform_q
        <
            group_degree_helper<Summand>,
            Groups
        >;
};

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
    bool Translation,
    typename Groups
>
struct error_expression
{
private:
    static constexpr Real delta_1 =
        static_filter_error<Expression, double>::error;
    using expanded = expand_polynomial<Expression, Translation>;
    using group_degrees = typename group_degrees_helper<Groups>::template fn
        <
            boost::mp11::mp_front<expanded>
        >;
    using group_expressions = boost::mp11::mp_transform_q
        <
            group_expression_helper<expanded>,
            Groups
        >;
    using group_powers = boost::mp11::mp_transform
        <
            power,
            group_expressions,
            group_degrees
        >;
    using scale = boost::mp11::mp_fold
        <
            boost::mp11::mp_pop_front<group_powers>,
            boost::mp11::mp_front<group_powers>,
            product
        >;
    static constexpr std::size_t total_degree =
        boost::mp11::mp_fold
            <
                group_degrees,
                boost::mp11::mp_int<0>,
                boost::mp11::mp_plus
            >::value;
    static constexpr Real delta_1_cor =
        delta_1 * round_up_1_n<Real>(total_degree);
    struct delta_constant : public static_constant_interface<Real>
    {
        static constexpr Real value = delta_1_cor;
        static constexpr bool non_negative = true;
    };
public:
    using type = product<delta_constant, scale>;
};

template
<
    typename Expression
>
struct trivial_groups_impl
{
private:
    using expanded = expand_polynomial<Expression, false>;
    using all_args = boost::mp11::mp_flatten<expanded, expand_product<>>;
public:
    using type = fpg_groups
        <
            boost::mp11::mp_rename
            <
                boost::mp11::mp_unique<all_args>,
                fpg_group
            >
        >;
};

template<typename Expression> using trivial_groups =
    typename trivial_groups_impl<Expression>::type;

template
<
    typename Expression,
    typename Groups
>
struct is_exact_cover_impl
{
private:
    using all_args = boost::mp11::mp_front<trivial_groups<Expression>>;
    using group_union = boost::mp11::mp_flatten<Groups, fpg_group<>>;
public:
    using type = boost::mp11::mp_and
        <
            boost::mp11::mp_is_set<group_union>,
            boost::mp11::mp_same
                <
                    boost::mp11::mp_size<all_args>,
                    boost::mp11::mp_size<group_union>
                >,
            boost::mp11::mp_empty
                <
                    boost::mp11::mp_set_difference
                    <
                        all_args,
                        boost::mp11::mp_unique<group_union>
                    >
                >
        >;
};

template <typename Expression, typename Groups>
using is_exact_cover = typename is_exact_cover_impl<Expression, Groups>::type;

template <typename Expression, typename Groups>
struct translations_in_groups_impl
{
private:
    using expanded = expand_polynomial<Expression, true>;
    using all_leaves = boost::mp11::mp_flatten<expanded, expand_product<>>;
    using all_differences = boost::mp11::mp_unique
        <
            boost::mp11::mp_copy_if
            <
                all_leaves,
                is_difference
            >
        >;
    template <typename Diff>
    struct group_contains_helper
    {
        template <typename Group>
        using fn = is_in_group<Diff, Group>;
    };
    template <typename Diff>
    using groups_contain = boost::mp11::mp_any_of_q
        <
            Groups,
            group_contains_helper<Diff>
        >;
public:
    using type = boost::mp11::mp_all_of
        <
            all_differences,
            groups_contain
        >;
};

template <typename Expression, bool Translation, typename Groups>
struct group_degree_equal_in_all_summands_impl
{
private:
    using expanded = expand_polynomial<Expression, Translation>;
    using all_group_degrees = boost::mp11::mp_transform_q
        <
            group_degrees_helper<Groups>,
            expanded
        >;
public:
    using type = boost::mp11::mp_bool
        <
            boost::mp11::mp_size
                <
                    boost::mp11::mp_unique<all_group_degrees>
                >::value == 1
        >;
};

template <typename Expression, bool Translation, typename Groups>
using valid_groups = boost::mp11::mp_and
    <
        typename group_degree_equal_in_all_summands_impl
            <
                Expression,
                Translation,
                Groups
            >::type,
        typename translations_in_groups_impl<Expression, Groups>::type,
        typename is_exact_cover_impl<Expression, Groups>::type
    >;

template <typename Expression, bool Translation>
struct auto_groups_impl
{
    using type = trivial_groups<Expression>;
};

template <typename Expression>
struct auto_groups_impl<Expression, true>
{
private:
    using translation_groups =
        typename translation_auto_groups<Expression>::type;
public:
    using type = boost::mp11::mp_if
        <
            valid_groups<Expression, true, translation_groups>,
            translation_groups,
            trivial_groups<Expression>
        >;
};

template <typename Expression, bool Translation>
using auto_groups = typename auto_groups_impl<Expression, Translation>::type;

} // fpg

template
<
    typename Expression,
    typename CalculationType,
    bool Translation = true,
    typename Groups = fpg::auto_groups<Expression, Translation>
>
using fpg_error_expression = typename fpg::error_expression
        <
            Expression,
            CalculationType,
            Translation,
            Groups
        >::type;

template
<
    typename Expression,
    typename CalculationType,
    bool Translation = true,
    typename Groups = fpg::auto_groups<Expression, Translation>
>
using fpg_semi_static = semi_static_filter
        <
            Expression,
            CalculationType,
            fpg_error_expression
                <
                    Expression,
                    CalculationType,
                    Translation,
                    Groups
                >
        >;

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_FPG_ERROR_BOUND_HPP
