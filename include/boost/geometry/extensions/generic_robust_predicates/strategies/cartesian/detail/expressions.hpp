// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPRESSIONS_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPRESSIONS_HPP

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template
<
    typename A11, typename A12,
    typename A21, typename A22
>
using det2x2 = difference
    <
        product<A11, A22>,
        product<A12, A21>
    >;

using orient2d = det2x2
    <
        difference <_1, _5>, difference<_2, _6>,
        difference <_3, _5>, difference<_4, _6>
    >;

template
<
    typename A11, typename A12, typename A13,
    typename A21, typename A22, typename A23,
    typename A31, typename A32, typename A33
>
struct det3x3_impl
{
private:
    using minor1 = product<A11, det2x2<A22, A23, A32, A33>>;
    using minor2 = product<A21, det2x2<A12, A13, A32, A33>>;
    using minor3 = product<A31, det2x2<A12, A13, A22, A23>>;
public:
    using type = sum<difference<minor1, minor2>, minor3>;
};

template
<
    typename A11, typename A12, typename A13,
    typename A21, typename A22, typename A23,
    typename A31, typename A32, typename A33
>
using det3x3 = typename det3x3_impl
    <
        A11, A12, A13,
        A21, A22, A23,
        A31, A32, A33
    >::type;

using orient3d = det3x3
    <
        difference<_1, _10>, difference<_2, _11>, difference<_3, _12>,
        difference<_4, _10>, difference<_5, _11>, difference<_6, _12>,
        difference<_7, _10>, difference<_8, _11>, difference<_9, _12>
    >;

struct incircle_impl
{
private:
    using adx = difference<_1, _7>;
    using ady = difference<_2, _8>;
    using bdx = difference<_3, _7>;
    using bdy = difference<_4, _8>;
    using cdx = difference<_5, _7>;
    using cdy = difference<_6, _8>;
    using alift = sum<product<adx, adx>, product<ady, ady>>;
    using blift = sum<product<bdx, bdx>, product<bdy, bdy>>;
    using clift = sum<product<cdx, cdx>, product<cdy, cdy>>;
public:
    using type = det3x3
        <
            alift, adx, ady,
            blift, bdx, bdy,
            clift, cdx, cdy
        >;
};

using incircle = incircle_impl::type;

template
<
    typename A11, typename A12, typename A13, typename A14,
    typename A21, typename A22, typename A23, typename A24,
    typename A31, typename A32, typename A33, typename A34,
    typename A41, typename A42, typename A43, typename A44
>
struct det4x4_impl
{
private:
    using minor1 = product
        <
            A11,
            det3x3<A22, A23, A24, A32, A33, A34, A42, A43, A44>
        >;
    using minor2 = product
        <
            A21,
            det3x3<A12, A13, A14, A32, A33, A34, A42, A43, A44>
        >;
    using minor3 = product
        <
            A31,
            det3x3<A12, A13, A14, A22, A23, A24, A42, A43, A44>
        >;
    using minor4 = product
        <
            A41,
            det3x3<A12, A13, A14, A22, A23, A24, A32, A33, A34>
        >;
public:
    using type = sum
        <
            difference<minor1, minor2>,
            difference<minor3, minor4>
        >;
};

template
<
    typename A11, typename A12, typename A13, typename A14,
    typename A21, typename A22, typename A23, typename A24,
    typename A31, typename A32, typename A33, typename A34,
    typename A41, typename A42, typename A43, typename A44
>
using det4x4 = typename det4x4_impl
    <
        A11, A12, A13, A14,
        A21, A22, A23, A24,
        A31, A32, A33, A34,
        A41, A42, A43, A44
    >::type;

struct insphere_impl
{
private:
    using aex = difference< _1, _13>;
    using aey = difference< _2, _14>;
    using aez = difference< _3, _15>;
    using bex = difference< _4, _13>;
    using bey = difference< _5, _14>;
    using bez = difference< _6, _15>;
    using cex = difference< _7, _13>;
    using cey = difference< _8, _14>;
    using cez = difference< _9, _15>;
    using dex = difference<_10, _13>;
    using dey = difference<_11, _14>;
    using dez = difference<_12, _15>;
    using alift = sum
        <
            product<aex, aex>,
            sum<product<aey, aey>, product<aez, aez>>
        >;
    using blift = sum
        <
            product<bex, bex>,
            sum<product<bey, bey>, product<bez, bez>>
        >;
    using clift = sum
        <
            product<cex, cex>,
            sum<product<cey, cey>, product<cez, cez>>
        >;
    using dlift = sum
        <
            product<dex, dex>,
            sum<product<dey, dey>, product<dez, dez>>
        >;
public:
    using type = det4x4
        <
            aex, aey, aez, alift,
            bex, bey, bez, blift,
            cex, cey, cez, clift,
            dex, dey, dez, dlift
        >;
};

using insphere = insphere_impl::type;

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPRESSIONS_HPP
