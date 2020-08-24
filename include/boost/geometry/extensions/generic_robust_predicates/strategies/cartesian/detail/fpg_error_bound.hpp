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

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

namespace fpg
{

// The following is an overestimation of ulp by at most a factor of 2.
// This could be improved.
template <typename Real>
constexpr Real ulp(Real d)
{
    assert( d >= 0 );
    return d * std::numeric_limits<Real>::epsilon();
}

// The following corrects for not being able to rounding towards
// infinity at compile-time.
template <typename Real>
constexpr Real round_up(Real d)
{
    return d + ulp(d);
}

template
<
    typename Expression,
    typename Real,
    operator_types Op = Expression::operator_type,
    bool IsLeaf = Expression::is_leaf
>
struct static_filter_error
{
};

template <typename Expression, typename Real>
struct static_filter_error <Expression, Real, operator_types::no_op, true>
{
    static constexpr Real bound = 1.0;
    static constexpr Real error = 0.0;
};

template <typename Expression, typename Real>
struct static_filter_error <Expression, Real, operator_types::sum, false>
{
private:
    using left_sfe = static_filter_error<typename Expression::left, Real>;
    using right_sfe =
        static_filter_error<typename Expression::right, Real>;
    static constexpr Real bound_sum =
        round_up(left_sfe::bound + right_sfe::bound);
    static constexpr Real u = ulp<Real>(bound_sum) / 2;
public:
    static constexpr Real bound = round_up(bound_sum + u);
    static constexpr Real error =
        round_up(round_up(u + left_sfe::error) + right_sfe::error);
};

template <typename Expression, typename Real>
struct static_filter_error<Expression, Real, operator_types::difference, false>
{
private:
    using left_sfe = static_filter_error<typename Expression::left, Real>;
    using right_sfe =
        static_filter_error<typename Expression::right, Real>;
    static constexpr Real bound_sum =
        round_up(left_sfe::bound + right_sfe::bound);
    static constexpr Real u = ulp<Real>(bound_sum) / 2;
    static constexpr bool input_translation =
        Expression::left::is_leaf && Expression::right::is_leaf;
public:
    static constexpr Real bound = 
        input_translation ?
              1
            : round_up(bound_sum + u);
    static constexpr Real error =
        input_translation ?
              ulp(Real(1)) / 2
            : round_up(round_up(u + left_sfe::error) + right_sfe::error);
};

template <typename Expression, typename Real>
struct static_filter_error<Expression, Real, operator_types::product, false>
{
private:
    using left_sfe = static_filter_error<typename Expression::left, Real>;
    using right_sfe =
        static_filter_error<typename Expression::right, Real>;
    static constexpr Real bound_product = 
        round_up(left_sfe::bound * right_sfe::bound);
    static constexpr Real u = ulp<Real>(bound_product) / 2;
public:
    static constexpr Real bound = round_up( bound_product + u );
private:
    static constexpr Real ee = round_up(left_sfe::error * right_sfe::error);
    static constexpr Real eb = round_up(left_sfe::error * right_sfe::bound);
    static constexpr Real be = round_up(left_sfe::bound * right_sfe::error);
public:
    static constexpr Real error = round_up(  round_up(u + ee)
                                           + round_up(eb + be));
};

template
<
    typename Expression,
    bool Translation,
    operator_types Op = Expression::operator_type
>
struct decomposition_anchor_impl
{
    using type = boost::mp11::mp_bool<Expression::is_leaf>;
};

template
<
    typename Expression
>
struct decomposition_anchor_impl
    <
        Expression,
        true,
        operator_types::difference
    >
{
    using type = boost::mp11::mp_bool
        <
            Expression::left::is_leaf && Expression::right::is_leaf
        >;
};

template <typename Expression, bool Translation>
using decomposition_anchor =
    typename decomposition_anchor_impl<Expression, Translation>::type;

template
<
    typename Expression,
    bool Translation = true,
    operator_types Op = Expression::operator_type,
    typename Anchor = decomposition_anchor<Expression, Translation>
>
struct decompose_polynomial_impl {};

template <typename ...> struct decomp_sum {};
template <typename ...> struct decomp_product {};

template
<
    typename Expression,
    bool Translation,
    operator_types Op
>
struct decompose_polynomial_impl
    <
        Expression,
        Translation,
        Op,
        boost::mp11::mp_true
    >
{
    using type = decomp_sum<decomp_product<Expression>>;
};

template
<
    typename Expression,
    bool Translation
>
struct decompose_polynomial_impl
    <
        Expression,
        Translation,
        operator_types::product,
        boost::mp11::mp_false
    >
{
private:
    using ld = typename decompose_polynomial_impl
        <
            typename Expression::left,
            Translation
        >::type;
    using rd = typename decompose_polynomial_impl
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
struct decompose_polynomial_impl
    <
        Expression,
        Translation,
        operator_types::sum,
        boost::mp11::mp_false
    >
{
private:
    using ld = typename decompose_polynomial_impl
        <
            typename Expression::left,
            Translation
        >::type;
    using rd = typename decompose_polynomial_impl
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
struct decompose_polynomial_impl
    <
        Expression,
        Translation,
        operator_types::difference,
        boost::mp11::mp_false
    >
{
private:
    using ld = typename decompose_polynomial_impl
        <
            typename Expression::left,
            Translation
        >::type;
    using rd = typename decompose_polynomial_impl
        <
            typename Expression::right,
            Translation
        >::type;
public:
    using type = boost::mp11::mp_append<ld, rd>;
};

template <typename Expression, bool Translation = true>
using decompose_polynomial =
    typename decompose_polynomial_impl<Expression, Translation>::type;

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
    static_assert(left_degree::value == right_degree::value,
                  "Must not mix degrees in difference.");
public:
    using type = left_degree;
};

template <typename LeafOrDiff, typename DegreeMap>
using degree = typename degree_impl<LeafOrDiff, DegreeMap>::type;

template <typename DecompProduct, typename DegreeMap>
struct product_degree_impl
{
private:
    using degree_q = boost::mp11::mp_bind_back<degree, DegreeMap>;
    using degrees = boost::mp11::mp_transform_q<degree_q, DecompProduct>;
public:
    using type = boost::mp11::mp_fold
        <
            degrees,
            boost::mp11::mp_int<0>,
            boost::mp11::mp_plus
        >;
};

template <typename DecompProduct, typename DegreeMap>
using product_degree =
    typename product_degree_impl<DecompProduct, DegreeMap>::type;

template <typename DecompPolynomial, typename DegreeMap>
struct decomp_polynomial_degree_impl
{
public:
    using type = typename product_degree_impl
        <
            boost::mp11::mp_first<DecompPolynomial>,
            DegreeMap
        >::type;
private:
    using degree_q = boost::mp11::mp_bind_back<product_degree, DegreeMap>;
    using degrees = boost::mp11::mp_transform_q<degree_q, DecompPolynomial>;
    template <typename Degree>
    using correct = boost::mp11::mp_bool<Degree::value == type::value>;
    using all_correct = boost::mp11::mp_all_of<degrees, correct>;
    static_assert(all_correct::value, "Polynomial must be homogenous.");
};

template <typename DecompPolynomial>
struct translation_auto_group
{
private:
    template <typename T>
    using is_difference =
        boost::mp_bool<T::operator_type == operator_typs::difference>;
    using translations =
        boost::mp11::mp_unique
            <
                boost::mp11::mp_copy_if
                    <
                        boost::mp11::mp_flatten
                            <
                                DecompPolynomial,
                                decomp_product<>
                            >,
                        is_difference
                    >
            >;
    //TODO: to be continued...
};

} // fpg

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_FPG_ERROR_BOUND_HPP
