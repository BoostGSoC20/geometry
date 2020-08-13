// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPANSION_ARITHMETIC_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPANSION_ARITHMETIC_HPP

#include <cassert>
#include <cmath>
#include <array>
#include <algorithm>
#include <iterator>
#include <type_traits>

#include <boost/mp11/integral.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/utility.hpp>
#include <boost/mp11/algorithm.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

//TODO: Reevaluate the thresholds for various summation algorithms
//TODO: Make use of zero-elimination
//TODO: Reevaluate the zero-elimination threshold
//TODO: Evaluate the use of AVX2-Gather instructions for post-summation zero-elimination

struct abs_comp
{
    template <typename Real>
    constexpr bool operator()(Real a, Real b) const
    {
        return std::abs(a) < std::abs(b);
    }
};

template <bool Negate, typename Real>
struct negate_impl
{
    static constexpr Real apply(Real a)
    {
        return a;
    }
};

template <typename Real>
struct negate_impl<true, Real>
{
    static constexpr Real apply(Real a)
    {
        return -a;
    }
};

template <bool Negate, typename Real>
constexpr Real negate(Real a)
{
    return negate_impl<Negate, Real>::apply(a);
}

namespace debug_expansion
{

unsigned long long round_to_power_of_two(unsigned long long num)
{
    if (!num)
    {
        return 0;
    }
    unsigned long long ret = 1;
    while (num >>= 1)
    {
        ret <<= 1;
    }
    return ret;
}

template <typename Real>
inline bool nonoverlapping(Real a, Real b)
{
    int a_exp, b_exp, min_exp, max_exp;
    Real a_mant = std::frexp(a, &a_exp);
    Real b_mant = std::frexp(b, &b_exp);
    unsigned long long scale = 1ULL<<63;
    unsigned long long a_mantll = std::abs(a_mant) * scale;
    unsigned long long b_mantll = std::abs(b_mant) * scale;
    if(a_mantll == 0 || b_mantll == 0)
    {
        return true;
    }
    unsigned long long min_mantll, max_mantll;
    if(a_exp < b_exp)
    {
        min_exp = a_exp;
        max_exp = b_exp;
        min_mantll = a_mantll;
        max_mantll = b_mantll;
    }
    else
    {
        min_exp = b_exp;
        max_exp = a_exp;
        min_mantll = b_mantll;
        max_mantll = a_mantll;
    }
    int scale_down = max_exp - min_exp;
    if(scale_down > std::numeric_limits<Real>::digits)
    {
        return true;
    }
    unsigned long long min_mantll_sc = min_mantll >> scale_down;
    auto min_mantll_sc_rd = round_to_power_of_two(min_mantll_sc);
    bool result = (max_mantll % ( 2 * min_mantll_sc_rd )) == 0;

    return result;
}

template <typename Real>
inline bool nonadjacent(Real a, Real b)
{
    bool t1 = nonoverlapping(a, b);
    bool t2 = nonoverlapping(a, 2*b);
    bool t3 = nonoverlapping(2*a, b);
    return t1 && t2 && t3;
}

template <typename Iter>
inline bool expansion_nonoverlapping(Iter begin, Iter end)
{
    auto lesser = *begin;
    for(auto it = begin + 1; it < end; ++it)
    {
        if( *it != 0)
        {
            if(lesser > std::abs(*it) || !nonoverlapping(lesser, *it) )
            {
                return false;
            }
            lesser = *it;
        }
    }
    return true;
}

template <typename Iter>
inline bool expansion_nonadjacent(Iter begin, Iter end)
{
    auto lesser = *begin;
    for(auto it = begin + 1; it < end; ++it)
    {
        if( *it != 0)
        {
            if(lesser > std::abs(*it) || !nonadjacent(lesser, *it))
            {
                return false;
            }
            lesser = *it;
        }
    }
    return true;
}

template <typename Iter>
inline bool expansion_strongly_nonoverlapping(Iter begin, Iter end)
{
    using Real = typename std::iterator_traits<Iter>::value_type;
    Real lesser = *begin;
    Real previous = 0.0;
    for(auto it = begin + 1; it < end; ++it)
    {
        if( *it != 0)
        {
            if(lesser > std::abs(*it) || !nonoverlapping(lesser, *it))
            {
                return false;
            }
            if( !nonadjacent(lesser, *it) )
            {
                int exp_now, exp_lesser;
                double lesser_mant = std::frexp(lesser, &exp_lesser);
                double now_mant = std::frexp(*it, &exp_now);
                if(    std::abs(std::frexp(lesser, &exp_lesser)) != 0.5
                    || std::abs(std::frexp(*it, &exp_now)) != 0.5 )
                {
                    return false;
                }
                if( !nonadjacent( lesser, previous ) )
                {
                    return false;
                }
            }
            previous = lesser;
            lesser = *it;
        }
    }
    return true;
}

} // namespace debug_expansion

template <typename Real>
constexpr Real two_sum_tail(Real a, Real b, Real x)
{
    Real b_virtual = x - a;
    Real a_virtual = x - b_virtual;
    Real b_rounded = b - b_virtual;
    Real a_rounded = a - a_virtual;
    Real y = a_rounded + b_rounded;

    assert(debug_expansion::nonadjacent(x, y));
    return y;
}

template <typename Real>
constexpr Real fast_two_sum_tail(Real a, Real b, Real x)
{
    assert(std::abs(a) >= std::abs(b) || a == 0);
    Real b_virtual = x - a;
    Real y = b - b_virtual;
    assert(debug_expansion::nonadjacent(x, y));
    return y;
}

template <typename Real>
constexpr Real two_difference_tail(Real a, Real b, Real x)
{
    Real b_virtual = a - x;
    Real a_virtual = x + b_virtual;
    Real b_rounded = b_virtual - b;
    Real a_rounded = a - a_virtual;
    Real y = a_rounded + b_rounded;
    assert(debug_expansion::nonadjacent(x, y));
    return y;
}

template <typename Real>
constexpr Real fast_two_difference_tail(Real a, Real b, Real x)
{
    assert(std::abs(a) >= std::abs(b) || a == 0);
    Real b_virtual = a - x;
    Real y = b_virtual - b;
    assert(debug_expansion::nonadjacent(x, y));
    return y;
}

template <typename Real>
constexpr Real two_product_tail(Real a, Real b, Real x)
{
    Real y = std::fma(a, b, -x);
    assert(debug_expansion::nonadjacent(x, y));
    return y;
}

template
<
    typename InIter,
    typename OutIter,
    typename Real,
    bool NegateE = false,
    bool NegateB = false
>
inline OutIter grow_expansion(InIter e_begin,
                              InIter e_end,
                              Real b,
                              OutIter h_begin,
                              OutIter h_end)
{
    assert(std::distance(e_begin, e_end) + 1 == std::distance(h_begin, h_end));
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));

    Real Q = negate<NegateB>(b);
    auto h_it = h_begin;
    for(auto e_it = e_begin; e_it != e_end; ++e_it)
    {
        Real Q_new = negate<NegateE>(*e_it) + Q;
        *h_it = two_sum_tail(negate<NegateE>(*e_it), Q, Q_new);
        Q = Q_new;
        ++h_it;
    }
    *h_it = Q;
    ++h_it;
    assert(debug_expansion::expansion_nonoverlapping(h_begin, h_end));
    assert(  !debug_expansion::expansion_nonadjacent(e_begin, e_end)
           || debug_expansion::expansion_nonadjacent(h_begin, h_end) );
    return h_it;
}

template
<
    typename InIter,
    typename OutIter,
    typename Real
>
inline OutIter grow_expansion_difference(InIter e_begin,
                                         InIter e_end,
                                         Real b,
                                         OutIter h_begin,
                                         OutIter h_end)
{
    return grow_expansion
        <
            InIter,
            OutIter,
            Real,
            false,
            true
        >(e_begin, e_end, b, h_begin, h_end);
}

template
<
    typename InIter,
    typename OutIter,
    typename Real,
    bool NegateE = false,
    bool NegateB = false
>
inline OutIter grow_expansion_ze(InIter e_begin,
                                 InIter e_end,
                                 Real b,
                                 OutIter h_begin,
                                 OutIter h_end)
{
    assert(std::distance(e_begin, e_end) + 1 <= std::distance(h_begin, h_end));
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));
    Real Q = negate<NegateB>(b);
    auto h_it = h_begin;
    for(auto e_it = e_begin; e_it != e_end; ++e_it)
    {
        Real Q_new = negate<NegateE>(*e_it) + Q;
        Real h_new = two_sum_tail(negate<NegateE>(*e_it), Q, Q_new);
        Q = Q_new;
        if(h_new != 0)
        {
            *h_it = h_new;
            ++h_it;
        }
    }
    if( Q != 0 || h_it == h_begin )
    {
        *h_it = Q;
        ++h_it;
    }
    assert(debug_expansion::expansion_nonoverlapping(h_begin, h_it));
    assert(  !debug_expansion::expansion_nonadjacent(e_begin, e_end)
           || debug_expansion::expansion_nonadjacent(h_begin, h_it));
    return h_it;
}

template
<
    typename InIter,
    typename OutIter,
    typename Real
>
inline OutIter grow_expansion_difference_ze(InIter e_begin,
                                            InIter e_end,
                                            Real b,
                                            OutIter h_begin,
                                            OutIter h_end)
{
    return grow_expansion_ze
        <
            InIter,
            OutIter,
            Real,
            false,
            true
        >(e_begin, e_end, b, h_begin, h_end);
}

template
<
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool NegateE = false,
    bool NegateF = false
>
inline OutIter expansion_sum(InIter1 e_begin,
                             InIter1 e_end,
                             InIter2 f_begin,
                             InIter2 f_end,
                             OutIter h_begin,
                             OutIter h_end)
{
    assert(std::distance(e_begin, e_end) + std::distance(f_begin, f_end)
            == std::distance(h_begin, h_end));
    using Real = typename std::iterator_traits<InIter1>::value_type;
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));
    assert(debug_expansion::expansion_nonoverlapping(f_begin, f_end));
    const auto elen = std::distance(e_begin, e_end);
    const auto flen = std::distance(f_begin, f_end);
    auto f_it = f_begin;
    grow_expansion
        <
            InIter1,
            OutIter,
            Real,
            NegateE,
            NegateF
        >(e_begin, e_end, *f_it, h_begin, h_begin + elen + 1);
    ++f_it;
    for(auto i = 1; i < flen; ++i)
    {
        grow_expansion
            <
                InIter1,
                OutIter,
                Real,
                false,
                NegateF
            >(h_begin + i,
              h_begin + i + elen,
              *f_it,
              h_begin + i,
              h_begin + i + elen + 1);
        ++f_it;
    }
    assert(debug_expansion::expansion_nonoverlapping(h_begin, h_end));
    assert(   !debug_expansion::expansion_nonadjacent(e_begin, e_end)
           || !debug_expansion::expansion_nonadjacent(f_begin, f_end)
           || debug_expansion::expansion_nonadjacent(h_begin, h_begin + elen + flen));
    return h_begin + elen + flen;
}

template
<
    typename InIter1,
    typename InIter2,
    typename OutIter
>
inline OutIter expansion_difference(InIter1 e_begin,
                                    InIter1 e_end,
                                    InIter2 f_begin,
                                    InIter2 f_end,
                                    OutIter h_begin,
                                    OutIter h_end)
{
    return expansion_sum
        <
            InIter1,
            InIter2,
            OutIter,
            false,
            true
        >(e_begin, e_end, f_begin, f_end, h_begin, h_end);
}

template
<
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool NegateE = false,
    bool NegateF = false
>
inline OutIter expansion_sum_ze(InIter1 e_begin,
                                InIter1 e_end,
                                InIter2 f_begin,
                                InIter2 f_end,
                                OutIter h_begin,
                                OutIter h_end)
{
    assert(std::distance(e_begin, e_end) + std::distance(f_begin, f_end)
            == std::distance(e_begin, e_end));
    using Real = typename std::iterator_traits<InIter1>::value_type;
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));
    assert(debug_expansion::expansion_nonoverlapping(f_begin, f_end));
    const auto elen = std::distance(e_begin, e_end);
    const auto flen = std::distance(f_begin, f_end);
    auto f_it = f_begin;
    auto h_end_new = grow_expansion_ze
        <
            InIter1,
            OutIter,
            Real,
            NegateE,
            NegateF
        >(e_begin, e_end, *f_it, h_begin, h_begin + elen + 1);
    ++f_it;
    for(auto i = 1; i < flen; ++i)
    {
        h_end_new = grow_expansion_ze
            <
                InIter1,
                OutIter,
                Real,
                false,
                NegateF
            >(h_begin + i,
              h_end_new,
              *f_it,
              h_begin + i,
              h_begin + h_end_new + 1);
        ++f_it;
    }
    assert(debug_expansion::expansion_nonoverlapping(h_begin, h_end_new));
    assert(   !debug_expansion::expansion_nonadjacent(e_begin, e_end)
           || !debug_expansion::expansion_nonadjacent(f_begin, f_end)
           || debug_expansion::expansion_nonadjacent(h_begin, h_end_new));
    return h_end_new;
}

template
<
    typename InIter1,
    typename InIter2,
    typename OutIter
>
inline OutIter expansion_difference_ze(InIter1 e_begin,
                                       InIter1 e_end,
                                       InIter2 f_begin,
                                       InIter2 f_end,
                                       OutIter h_begin,
                                       OutIter h_end)
{
    return expansion_sum_ze
        <
            InIter1,
            InIter2,
            OutIter,
            false,
            true
        >(e_begin, e_end, f_begin, f_end, h_begin, h_end);
}

template
<
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool NegateE = false,
    bool NegateF = false
>
inline OutIter fast_expansion_sum_inplace(InIter1 e_begin,
                                          InIter1 e_end,
                                          InIter2 f_begin,
                                          InIter2 f_end,
                                          OutIter h_begin,
                                          OutIter
                                          h_end)
{
    assert(e_end == f_begin);
    assert(h_begin == e_begin && h_end == f_end);
    assert(f_begin != h_begin);
    assert(std::distance(e_begin, e_end) + std::distance(f_begin, f_end)
            == std::distance(h_begin, h_end));
    using Real = typename std::iterator_traits<InIter1>::value_type;
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));
    assert(debug_expansion::expansion_nonoverlapping(f_begin, f_end));
    if(NegateE)
    {
        for(auto e_it = e_begin; e_it != e_end; ++e_it)
        {
            *e_it = -*e_it;
        }
    }
    if(NegateF)
    {
        for(auto f_it = f_begin; f_it != f_end; ++f_it)
        {
            *f_it = -*f_it;
        }
    }
    //TODO: The following may seem wasteful but in the more relevant case of zero-elimination
    //      it should reduce to a single rotate.

    auto e_old_end = e_end;
    e_end = std::remove(e_begin, e_end, Real(0));
    auto e_shortened = e_old_end - e_end;
    std::rotate(e_end, f_begin, f_end);
    f_begin = e_end;
    auto f_old_end = f_end;
    f_end = f_end - e_shortened;
    f_end = std::remove(e_end, f_end, Real(0));
    std::fill(f_end, f_old_end, Real(0));

    std::inplace_merge(e_begin, e_end, f_end, abs_comp{});

    auto g_it = e_begin;
    auto g_end = f_end;
    auto h_it = h_begin;
    Real Q = *g_it + *(g_it + 1);
    *h_it = fast_two_sum_tail(*(g_it + 1), *(g_it), Q);
    g_it += 2;
    ++h_it;
    for(; g_it < g_end; ++g_it)
    {
        Real Q_new = Q + *g_it;
        *h_it = two_sum_tail(Q, *g_it, Q_new);
        Q = Q_new;
        ++h_it;
    }
    *h_it = Q;
    ++h_it;
    assert(debug_expansion::expansion_strongly_nonoverlapping(e_begin, f_end));
    return h_end;
}

template
<
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool NegateE = false,
    bool NegateF = false
>
inline OutIter fast_expansion_sum_inplace_ze(InIter1 e_begin,
                                             InIter1 e_end,
                                             InIter2 f_begin,
                                             InIter2 f_end,
                                             OutIter h_begin,
                                             OutIter h_end)
{
    assert(e_end == f_begin);
    assert(h_begin == e_begin && h_end == f_end);
    assert(f_begin != h_begin);
    assert(std::distance(e_begin, e_end) + std::distance(f_begin, f_end)
            <= std::distance(h_begin, h_end));
    using Real = typename std::iterator_traits<InIter1>::value_type;
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));
    assert(debug_expansion::expansion_nonoverlapping(f_begin, f_end));
    if(NegateE)
    {
        for(auto e_it = e_begin; e_it != e_end; ++e_it)
        {
            *e_it = -*e_it;
        }
    }
    if(NegateF)
    {
        for(auto f_it = f_begin; f_it != f_end; ++f_it)
        {
            *f_it = -*f_it;
        }
    }
    std::inplace_merge(e_begin, e_end, f_end, abs_comp{});
    auto g_it = e_begin;
    auto g_end = f_end;
    auto h_it = h_begin;
    Real Q = *g_it + *(g_it + 1);
    Real h_new = fast_two_sum_tail(*(g_it + 1), *(g_it), Q);
    g_it += 2;
    if(h_new != 0)
    {
        *h_it = h_new;
        ++h_it;
    }
    for(; g_it != g_end; ++g_it)
    {
        Real Q_new = Q + *g_it;
        h_new = two_sum_tail(Q, *g_it, Q_new);
        if(h_new != 0)
        {
            *h_it = h_new;
            ++h_it;
        }
        Q = Q_new;
    }
    if( Q != 0 || h_it == h_begin )
    {
        *h_it = Q;
        ++h_it;
    }
    assert(debug_expansion::expansion_strongly_nonoverlapping(e_begin, h_it));
    return h_it;
}

template
<
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool NegateE = false,
    bool NegateF = false
>
inline OutIter fast_expansion_sum_not_inplace(InIter1 e_begin,
                                              InIter1 e_end,
                                              InIter2 f_begin,
                                              InIter2 f_end,
                                              OutIter h_begin,
                                              OutIter h_end)
{
    assert(e_begin != h_begin);
    assert(f_begin != h_begin);
    assert(std::distance(e_begin, e_end) + std::distance(f_begin, f_end)
            == std::distance(h_begin, h_end));
    using Real = typename std::iterator_traits<InIter1>::value_type;
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));
    assert(debug_expansion::expansion_nonoverlapping(f_begin, f_end));

    auto e_it = e_begin;
    auto f_it = f_begin;
    Real Q;
    if(std::abs(*f_it) > std::abs(*e_it))
    {
        Q = negate<NegateE>(*e_it);
        ++e_it;
    }
    else
    {
        Q = negate<NegateF>(*f_it);
        ++f_it;
    }
    auto h_it = h_begin;
    if ((e_it != e_end) && (f_it != f_end))
    {
        Real Q_new;
        if (std::abs(*f_it) > std::abs(*e_it))
        {
            Q_new = negate<NegateE>(*e_it) + Q;
            auto h_new = fast_two_sum_tail(negate<NegateE>(*e_it), Q, Q_new);
            *h_it = h_new;
            ++e_it;
        }
        else
        {
            Q_new = negate<NegateF>(*f_it) + Q;
            auto h_new = fast_two_sum_tail(negate<NegateF>(*f_it), Q, Q_new);
            *h_it = h_new;
            ++f_it;
        }
        Q = Q_new;
        ++h_it;
        while((e_it != e_end) && (f_it != f_end))
        {
            if (std::abs(*f_it) > std::abs(*e_it))
            {
                Q_new = negate<NegateE>(*e_it) + Q;
                *h_it = two_sum_tail(negate<NegateE>(*e_it), Q, Q_new);
                ++e_it;
            }
            else
            {
                Q_new = negate<NegateF>(*f_it) + Q;
                *h_it = two_sum_tail(negate<NegateF>(*f_it), Q, Q_new);
                ++f_it;
            }
            Q = Q_new;
            ++h_it;
        }
    }
    while(e_it != e_end)
    {
        Real Q_new = negate<NegateE>(*e_it) + Q;
        *h_it = two_sum_tail(negate<NegateE>(*e_it), Q, Q_new);
        Q = Q_new;
        ++e_it;
        ++h_it;
    }
    while(f_it != f_end)
    {
        Real Q_new = negate<NegateF>(*f_it) + Q;
        *h_it = two_sum_tail(negate<NegateF>(*f_it), Q, Q_new);
        Q = Q_new;
        ++f_it;
        ++h_it;
    }
    *h_it = Q;
    ++h_it;
    assert(debug_expansion::expansion_strongly_nonoverlapping(h_begin, h_end));
    return h_end;
}

template
<
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool inplace,
    bool NegateE,
    bool NegateF
>
struct fast_expansion_sum_impl
{
    static inline OutIter apply(InIter1 e_begin,
                                InIter1 e_end,
                                InIter2 f_begin,
                                InIter2 f_end,
                                OutIter h_begin,
                                OutIter h_end)
    {
        return fast_expansion_sum_not_inplace
            <
                InIter1,
                InIter2,
                OutIter,
                NegateE,
                NegateF
            >(e_begin, e_end, f_begin, f_end, h_begin, h_end);
    }
};

template
<
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool NegateE,
    bool NegateF
>
struct fast_expansion_sum_impl
    <
        InIter1,
        InIter2,
        OutIter,
        true,
        NegateE,
        NegateF
    >
{
    static inline OutIter apply(InIter1 e_begin,
                                InIter1 e_end,
                                InIter2 f_begin,
                                InIter2 f_end,
                                OutIter h_begin,
                                OutIter h_end)
    {
        return fast_expansion_sum_inplace
            <
                InIter1,
                InIter2,
                OutIter,
                NegateE,
                NegateF
            >(e_begin, e_end, f_begin, f_end, h_begin, h_end);
    }
};

template
<
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool inplace,
    bool NegateE,
    bool NegateF
>
inline OutIter fast_expansion_sum(InIter1 e_begin,
                                  InIter1 e_end,
                                  InIter2 f_begin,
                                  InIter2 f_end,
                                  OutIter h_begin,
                                  OutIter h_end)
{
    return fast_expansion_sum_impl
        <
            InIter1,
            InIter2,
            OutIter,
            inplace,
            NegateE,
            NegateF
        >::apply(e_begin, e_end, f_begin, f_end, h_begin, h_end);
}

template
<
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool inplace
>
inline OutIter fast_expansion_difference(InIter1 e_begin,
                                         InIter1 e_end,
                                         InIter2 f_begin,
                                         InIter2 f_end,
                                         OutIter h_begin,
                                         OutIter h_end)
{
    return fast_expansion_sum
        <
            InIter1,
            InIter2,
            OutIter,
            inplace,
            false,
            true
        >(e_begin, e_end, f_begin, f_end, h_begin, h_end);
}

template
<
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool NegateE = false,
    bool NegateF = false
>
inline OutIter fast_expansion_sum_not_inplace_ze(InIter1 e_begin,
                                                 InIter1 e_end,
                                                 InIter2 f_begin,
                                                 InIter2 f_end,
                                                 OutIter h_begin,
                                                 OutIter h_end)
{
    assert(std::distance(e_begin, e_end) + std::distance(f_begin, f_end)
            == std::distance(h_begin, h_end));
    using Real = typename std::iterator_traits<InIter1>::value_type;
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));
    assert(debug_expansion::expansion_nonoverlapping(f_begin, f_end));

    auto e_it = e_begin;
    auto f_it = f_begin;
    Real Q;
    if(std::abs(*f_it) > std::abs(*e_it))
    {
        Q = negate<NegateE>(*e_it);
        ++e_it;
    }
    else
    {
        Q = negate<NegateF>(*f_it);
        ++f_it;
    }
    auto h_it = h_begin;
    if ((e_it != e_end) && (f_it != f_end))
    {
        Real Q_new;
        Real h_new;
        if (std::abs(*f_it) > std::abs(*e_it))
        {
            Q_new = negate<NegateE>(*e_it) + Q;
            h_new = fast_two_sum_tail(negate<NegateE>(*e_it), Q, Q_new);
            ++e_it;
        }
        else
        {
            Q_new = negate<NegateF>(*f_it) + Q;
            h_new = fast_two_sum_tail(negate<NegateF>(*f_it), Q, Q_new);
            ++f_it;
        }
        Q = Q_new;
        if(h_new != 0)
        {
            *h_it = h_new;
            ++h_it;
        }
        while((e_it != e_end) && (f_it != f_end))
        {
            if (std::abs(*f_it) > std::abs(*e_it))
            {
                Q_new = negate<NegateE>(*e_it) + Q;
                h_new = two_sum_tail(negate<NegateE>(*e_it), Q, Q_new);
                ++e_it;
            }
            else
            {
                Q_new = negate<NegateF>(*f_it) + Q;
                h_new = two_sum_tail(negate<NegateF>(*f_it), Q, Q_new);
                ++f_it;
            }
            Q = Q_new;
            if(h_new != 0)
            {
                *h_it = h_new;
                ++h_it;
            }
        }
    }
    while(e_it != e_end)
    {
        Real Q_new = negate<NegateE>(*e_it) + Q;
        Real h_new = two_sum_tail(negate<NegateE>(*e_it), Q, Q_new);
        Q = Q_new;
        ++e_it;
        if(h_new != 0)
        {
            *h_it = h_new;
            ++h_it;
        }
    }
    while(f_it != f_end)
    {
        Real Q_new = negate<NegateF>(*f_it) + Q;
        Real h_new = two_sum_tail(negate<NegateF>(*f_it), Q, Q_new);
        Q = Q_new;
        ++f_it;
        if(h_new != 0)
        {
            *h_it = h_new;
            ++h_it;
        }
    }
    if( Q != 0 || h_it == h_begin)
    {
        *h_it = Q;
        ++h_it;
    }
    assert(debug_expansion::expansion_strongly_nonoverlapping(h_begin, h_it));
    return h_it;
}

template
<
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool inplace,
    bool NegateE,
    bool NegateF
>
struct fast_expansion_sum_ze_impl
{
    static inline OutIter apply(InIter1 e_begin,
                                InIter1 e_end,
                                InIter2 f_begin,
                                InIter2 f_end,
                                OutIter h_begin,
                                OutIter h_end)
    {
        return fast_expansion_sum_not_inplace_ze
            <
                InIter1,
                InIter2,
                OutIter,
                NegateE,
                NegateF
            >(e_begin, e_end, f_begin, f_end, h_begin, h_end);
    }
};

template
<
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool NegateE,
    bool NegateF
>
struct fast_expansion_sum_ze_impl
    <
        InIter1,
        InIter2,
        OutIter,
        true,
        NegateE,
        NegateF
    >
{
    static inline OutIter apply(InIter1 e_begin,
                                InIter1 e_end,
                                InIter2 f_begin,
                                InIter2 f_end,
                                OutIter h_begin,
                                OutIter h_end)
    {
        return fast_expansion_sum_inplace_ze
            <
                InIter1,
                InIter2,
                OutIter,
                NegateE,
                NegateF
            >(e_begin, e_end, f_begin, f_end, h_begin, h_end);
    }
};

template
<
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool inplace,
    bool NegateE,
    bool NegateF
>
inline OutIter fast_expansion_sum_ze(InIter1 e_begin,
                                     InIter1 e_end,
                                     InIter2 f_begin,
                                     InIter2 f_end,
                                     OutIter h_begin,
                                     OutIter h_end)
{
    return fast_expansion_sum_ze_impl
        <
            InIter1,
            InIter2,
            OutIter,
            inplace,
            NegateE,
            NegateF
        >::apply(e_begin, e_end, f_begin, f_end, h_begin, h_end);
}

template
<
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool inplace
>
inline OutIter fast_expansion_difference_ze(InIter1 e_begin,
                                            InIter1 e_end,
                                            InIter2 f_begin,
                                            InIter2 f_end,
                                            OutIter h_begin,
                                            OutIter h_end)
{
    return fast_expansion_sum_ze
        <
            InIter1,
            InIter2,
            OutIter,
            inplace,
            false,
            true
        >(e_begin, e_end, f_begin, f_end, h_begin, h_end);
}

template
<
    typename InIter,
    typename Real,
    typename OutIter
>
inline OutIter scale_expansion(InIter e_begin,
                               InIter e_end,
                               Real b,
                               OutIter h_begin,
                               OutIter h_end)
{
    assert(std::distance(e_begin, e_end) * 2 == std::distance(h_begin, h_end));
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));

    auto e_it = e_begin;
    auto h_it = h_begin;
    Real Q = *e_it * b;
    *h_it = two_product_tail(*e_it, b, Q);
    ++e_it;
    ++h_it;
    for(; e_it != e_end; ++e_it)
    {
        Real product_1 = *e_it * b;
        Real product_0 = two_product_tail(*e_it, b, product_1);
        Real sum = Q + product_0;
        *h_it = two_sum_tail(Q, product_0, sum);
        ++h_it;
        Q = product_1 + sum;
        *h_it = two_sum_tail(product_1, sum, Q);
        ++h_it;
    }
    *h_it = Q;
    ++h_it;

    assert( debug_expansion::expansion_nonoverlapping(h_begin, h_end) );
    assert(  !debug_expansion::expansion_nonadjacent(e_begin, e_end)
           || debug_expansion::expansion_nonadjacent(h_begin, h_end) );
    assert(  !debug_expansion::expansion_strongly_nonoverlapping(e_begin, e_end)
           || debug_expansion::expansion_strongly_nonoverlapping(h_begin, h_end) );
    return h_it;
}

template
<
    typename InIter,
    typename Real,
    typename OutIter
>
inline OutIter scale_expansion_ze(InIter e_begin,
                                  InIter e_end,
                                  Real b,
                                  OutIter h_begin,
                                  OutIter h_end)
{
    assert(std::distance(e_begin, e_end) * 2 == std::distance(h_begin, h_end));
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));
    auto e_it = e_begin;
    auto h_it = h_begin;
    Real Q = *e_it * b;
    Real h_new = two_product_tail(*e_it, b, Q);
    if(h_new != 0)
    {
        *h_it = h_new;
        ++h_it;
    }
    ++e_it;
    for(; e_it != e_end; ++e_it)
    {
        Real product_1 = *e_it * b;
        Real product_0 = two_product_tail(*e_it, b, product_1);
        Real sum = Q + product_0;
        Real h_new = two_sum_tail(Q, product_0, sum);
        if(h_new != 0)
        {
            *h_it = h_new;
            ++h_it;
        }
        Q = product_1 + sum;
        h_new = two_sum_tail(product_1, sum, Q);
        if(h_new != 0)
        {
            *h_it = h_new;
            ++h_it;
        }
    }
    if( Q != 0 || h_it == h_begin )
    {
        *h_it = Q;
        ++h_it;
    }
    assert(  !debug_expansion::expansion_nonadjacent(e_begin, e_end)
           || debug_expansion::expansion_nonadjacent(h_begin, h_it) );
    assert(  !debug_expansion::expansion_strongly_nonoverlapping(e_begin, e_end)
           || debug_expansion::expansion_strongly_nonoverlapping(h_begin, h_it) );
    assert( debug_expansion::expansion_nonoverlapping(h_begin, h_it) );
    return h_it;
}

template <typename Rhs>
struct greater_than_or_equal
{
    template<typename Lhs>
    using fn = boost::mp11::mp_bool< Lhs::value >= Rhs::value >;
};

template
<
    std::size_t e_length,
    std::size_t f_length,
    bool inplace = false,
    bool e_negate = false,
    bool f_negate = false,
    std::size_t h_length = e_length + f_length
>
struct expansion_plus_impl
{
    template<typename InIter1, typename InIter2, typename OutIter>
    static inline OutIter apply(
        InIter1 e_begin,
        InIter1 e_end,
        InIter2 f_begin,
        InIter2 f_end,
        OutIter h_begin,
        OutIter h_end)
    {
        return fast_expansion_sum
            <
                InIter1,
                InIter2,
                OutIter,
                inplace,
                e_negate,
                f_negate
            >(e_begin, e_end, f_begin, f_end, h_begin, h_end);
    }
};

template <bool inplace, bool e_negate, bool f_negate>
struct expansion_plus_impl<1, 1, inplace, e_negate, f_negate, 2>
{
    template<typename InIter1, typename InIter2, typename OutIter>
    static inline OutIter apply(InIter1 e_begin,
                                InIter1 e_end,
                                InIter2 f_begin,
                                InIter2 f_end,
                                OutIter h_begin,
                                OutIter h_end)
    {
        auto x = negate<e_negate>(*e_begin) + negate<f_negate>(*f_begin);
        auto y = two_sum_tail(negate<e_negate>(*e_begin),
                              negate<f_negate>(*f_begin), x);
        *h_begin = y;
        *(h_begin + 1) = x;
        return h_begin + 2;
    }
};

template
<
    std::size_t e_length,
    bool inplace,
    bool e_negate,
    bool f_negate,
    std::size_t h_length
>
struct expansion_plus_impl<e_length, 1, inplace, e_negate, f_negate, h_length>
{
    template<typename InIter1, typename InIter2, typename OutIter>
    static inline OutIter apply(InIter1 e_begin,
                                InIter1 e_end,
                                InIter2 f_begin,
                                InIter2 f_end,
                                OutIter h_begin,
                                OutIter h_end)
    {
        return grow_expansion
            <
                InIter1,
                OutIter,
                typename std::iterator_traits<InIter2>::value_type,
                e_negate,
                f_negate
            >(e_begin, e_end, *f_begin, h_begin, h_end);
    }
};

template
<
    std::size_t f_length,
    bool inplace,
    bool e_negate,
    bool f_negate,
    std::size_t h_length
>
struct expansion_plus_impl
    <
        1,
        f_length,
        inplace,
        e_negate,
        f_negate,
        h_length
    >
{
    template<typename InIter1, typename InIter2, typename OutIter>
    static inline OutIter apply(InIter1 e_begin,
                                InIter1 e_end,
                                InIter2 f_begin,
                                InIter2 f_end,
                                OutIter h_begin,
                                OutIter h_end)
    {
        return grow_expansion
            <
                InIter2,
                OutIter,
                typename std::iterator_traits<InIter1>::value_type,
                e_negate,
                f_negate
            >(f_begin, f_end, *e_begin, h_begin, h_end);
    }
};

template<bool inplace, bool e_negate, bool f_negate>
struct expansion_plus_impl<2, 2, inplace, e_negate, f_negate, 4>
{
    template<typename InIter1, typename InIter2, typename OutIter>
    static inline OutIter apply(InIter1 e_begin,
                                InIter1 e_end,
                                InIter2 f_begin,
                                InIter2 f_end,
                                OutIter h_begin,
                                OutIter h_end)
    {
        return expansion_sum
            <
                InIter1,
                InIter2,
                OutIter,
                e_negate,
                f_negate
            >(e_begin, e_end, f_begin, f_end, h_begin, h_end);
    }
};

template
<
    std::size_t e_length,
    std::size_t f_length,
    bool inplace,
    typename InIter1,
    typename InIter2,
    typename OutIter,
    std::size_t result = e_length + f_length
>
inline OutIter expansion_plus(InIter1 e_begin,
                              InIter1 e_end,
                              InIter2 f_begin,
                              InIter2 f_end,
                              OutIter h_begin,
                              OutIter h_end)
{
    return expansion_plus_impl<e_length, f_length, inplace>::
        apply(e_begin, e_end, f_begin, f_end, h_begin, h_end);
}

template
<
    std::size_t e_length,
    std::size_t f_length,
    bool inplace,
    typename InIter,
    typename Real,
    typename OutIter,
    std::size_t result = e_length + f_length
>
inline OutIter expansion_plus(InIter e_begin,
                              InIter e_end,
                              Real f,
                              OutIter h_begin,
                              OutIter h_end)
{
    static_assert( f_length == 1, "f_length must be 1 if f is a single component." );
    return grow_expansion(e_begin, e_end, f, h_begin, h_end);
}

template
<
    std::size_t e_length,
    std::size_t f_length,
    bool inplace,
    typename InIter,
    typename Real,
    typename OutIter,
    std::size_t result = e_length + f_length
>
inline OutIter expansion_plus(Real e,
                              InIter f_begin,
                              InIter f_end,
                              OutIter h_begin,
                              OutIter h_end)
{
    static_assert( e_length == 1, "e_length must be 1 if e is a single component." );
    return grow_expansion(f_begin, f_end, e, h_begin, h_end);
}

template
<
    std::size_t e_length,
    std::size_t f_length,
    bool inplace,
    typename Real,
    typename OutIter,
    std::size_t result = e_length + f_length
>
inline OutIter expansion_plus(
    Real e,
    Real f,
    OutIter h_begin,
    OutIter h_end)
{
    static_assert( f_length == 1 && e_length == 1, "e_length and f_length must be 1 if they are single components." );
    *(h_begin + 1) = e + f;
    *h_begin = two_sum_tail(e, f, *(h_begin + 1));
    return h_begin + 2;
}

template
<
    std::size_t e_length,
    std::size_t f_length,
    bool inplace,
    bool StageB,
    typename InIter1,
    typename InIter2,
    typename OutIter,
    std::size_t result = e_length + f_length
>
inline OutIter expansion_minus(InIter1 e_begin,
                               InIter1 e_end,
                               InIter2 f_begin,
                               InIter2 f_end,
                               OutIter h_begin,
                               OutIter h_end)
{
    return expansion_plus_impl<e_length, f_length, inplace, false, true>::
        apply(e_begin, e_end, f_begin, f_end, h_begin, h_end);
}

template
<
    std::size_t e_length,
    std::size_t f_length,
    bool inplace,
    bool StageB,
    typename InIter,
    typename Real,
    typename OutIter,
    std::size_t result = e_length + f_length
>
inline OutIter expansion_minus(InIter e_begin,
                               InIter e_end,
                               Real f,
                               OutIter h_begin,
                               OutIter h_end)
{
    static_assert(f_length == 1, "f_length must be 1 if f is a single component.");
    return expansion_plus
        <
            e_length,
            f_length,
            inplace
        >(e_begin, e_end, -f, h_begin, h_end);
}

template
<
    std::size_t e_length,
    std::size_t f_length,
    bool inplace,
    bool StageB,
    typename Real,
    typename InIter,
    typename OutIter,
    std::size_t result = e_length + f_length
>
inline OutIter expansion_minus(Real e,
                               InIter f_begin,
                               InIter f_end,
                               OutIter h_begin,
                               OutIter h_end)
{
    static_assert(e_length == 1, "e_length must be 1 if e is a single component.");
    return grow_expansion
        <
            InIter,
            OutIter,
            Real,
            false,
            true
        >(f_begin, f_end, e, h_begin, h_end);
}

template
<
    std::size_t e_length,
    std::size_t f_length,
    bool inplace,
    bool StageB,
    typename Real,
    typename OutIter,
    std::size_t result = e_length + f_length
>
inline std::enable_if_t<!StageB, OutIter> expansion_minus(Real e,
                                                         Real f,
                                                         OutIter h_begin,
                                                         OutIter h_end)
{
    static_assert(e_length == 1 && f_length == 1, "e_length and f_length must be 1 if they are single components.");
    *(h_begin + 1) = e - f;
    *h_begin = two_difference_tail(e, f, *(h_begin + 1));
    return h_begin + 2;
}

template
<
    std::size_t e_length,
    std::size_t f_length,
    bool inplace,
    bool StageB,
    typename Real,
    typename OutIter,
    std::size_t result = e_length + f_length
>
inline std::enable_if_t<StageB, OutIter> expansion_minus(Real e,
                                                         Real f,
                                                         OutIter h_begin,
                                                         OutIter h_end)
{
    static_assert(e_length == 1 && f_length == 1, "e_length and f_length must be 1 if they are single components.");
    *(h_begin) = e - f;
    return h_begin + 1;
}

template <std::size_t N>
struct minus_n_impl
{
    template<typename V>
    using fn = boost::mp11::mp_size_t<V::value - N>;
};

template
<
    typename IndexList,
    std::size_t index_length = boost::mp11::mp_size<IndexList>::value
>
struct distillation_impl
{
private:
    static constexpr std::size_t length =
        boost::mp11::mp_back<IndexList>::value;
    using first_greater =
        boost::mp11::mp_find_if_q
            <
                IndexList,
                greater_than_or_equal<boost::mp11::mp_size_t< (length + 1) / 2 >>
            >;
    using not_first =
        boost::mp11::mp_if
            <
                boost::mp11::mp_bool<(first_greater::value > 0)>,
                first_greater,
                boost::mp11::mp_size_t<first_greater::value + 1>
            >;
    using first_half = boost::mp11::mp_take<IndexList, not_first>;
    static constexpr std::size_t first_length =
        boost::mp11::mp_back<first_half>::value;
    using second_half = boost::mp11::mp_transform_q
        <
            minus_n_impl<first_length>,
            boost::mp11::mp_drop<IndexList, not_first>
        >;
public:
    template <typename Iter>
    static inline Iter apply(Iter begin, Iter end)
    {
        //The folliwng divide and conquer scheme is a primitive heuristic
        //and it would probably be a good idea to explore possible optimizations.
        //
        //TODO: Consider removing zeroes at this stage by writing output of second sub-
        //distillation to end of first sub-distillation.
        assert(std::distance(begin, end) <= length);
        auto first_end = distillation_impl<first_half>::
            apply(begin, begin + first_length);
        auto second_end = distillation_impl<second_half>::
            apply(begin + first_length, end);

        return expansion_plus
            <
                first_length,
                length - first_length,
                true
            >(begin,
              first_end,
              begin + first_length,
              second_end,
              begin,
              end);
    }
};

template<typename IndexList>
struct distillation_impl<IndexList, 1>
{
    template<typename Iter>
    static inline Iter apply(Iter begin, Iter end)
    {
        return end;
    }
};

template<typename IndexList, typename Iter>
inline Iter distillation(Iter begin, Iter end)
{
    return distillation_impl<IndexList>::apply(begin, end);
}

template
<
    std::size_t e_length,
    std::size_t f_length,
    std::size_t result = 2 * e_length * f_length,
    bool e_smaller = e_length <= f_length
>
struct expansion_times_impl
{
    template<typename InIter1, typename InIter2, typename OutIter>
    static inline OutIter apply(InIter1 e_begin,
                                InIter1 e_end,
                                InIter2 f_begin,
                                InIter2 f_end,
                                OutIter h_begin,
                                OutIter h_end)
    {
        assert(e_begin != h_begin && f_begin != h_begin);

        //TODO: Evaluate zero-elimination for very short expansions before multiplication.

        auto h_it = h_begin;
        for(auto e_it = e_begin; e_it != e_end; ++e_it)
        {
            scale_expansion(f_begin, f_end, *e_it, h_it, h_it + 2 * f_length);
            h_it += 2 * f_length;
        }
        using index_list = boost::mp11::mp_repeat_c
            <
                boost::mp11::mp_list<boost::mp11::mp_size_t<2 * f_length>>,
                e_length
            >;
        using ps_index_list = boost::mp11::mp_partial_sum
            <
                index_list,
                boost::mp11::mp_size_t<0>,
                boost::mp11::mp_plus
            >;

        h_it = distillation<ps_index_list>(h_begin, h_end);
        assert(debug_expansion::expansion_nonoverlapping(h_begin, h_it));
        return h_it;
    }
};

template
<
    std::size_t e_length,
    std::size_t f_length,
    std::size_t result
>
struct expansion_times_impl<e_length, f_length, result, false>
{
    template<typename InIter1, typename InIter2, typename OutIter>
    static inline OutIter apply(InIter1 e_begin,
                                InIter1 e_end,
                                InIter2 f_begin,
                                InIter2 f_end,
                                OutIter h_begin,
                                OutIter h_end)
    {
        return expansion_times_impl<f_length, e_length>
            ::apply(f_begin,
                    f_end,
                    e_begin,
                    e_end,
                    h_begin,
                    h_end);
    }
};

template<>
struct expansion_times_impl<1, 1, 2, true>
{
    template<typename InIter1, typename InIter2, typename OutIter>
    static inline OutIter apply(InIter1 e_begin,
                                InIter1 e_end,
                                InIter2 f_begin,
                                InIter2 f_end,
                                OutIter h_begin,
                                OutIter h_end)
    {
        auto x = *e_begin * *f_begin;
        auto y = two_product_tail(*e_begin, *f_begin, x);
        *h_begin = y;
        *(h_begin + 1) = x;
        return h_begin + 2;
    }
};

template
<
    std::size_t e_length,
    std::size_t f_length,
    bool inplace,
    typename InIter1,
    typename InIter2,
    typename OutIter,
    std::size_t result = 2 * e_length * f_length
>
inline OutIter expansion_times(InIter1 e_begin,
                               InIter1 e_end,
                               InIter2 f_begin,
                               InIter2 f_end,
                               OutIter h_begin,
                               OutIter h_end)
{
    return expansion_times_impl<e_length, f_length>::
        apply(e_begin, e_end, f_begin, f_end, h_begin, h_end);
}

template
<
    std::size_t e_length,
    std::size_t f_length,
    bool inplace,
    typename InIter,
    typename Real,
    typename OutIter,
    std::size_t result = 2 * e_length
>
inline OutIter expansion_times(InIter e_begin,
                               InIter e_end,
                               Real f,
                               OutIter h_begin,
                               OutIter h_end)
{
    static_assert(f_length == 1, "f_length must be 1 if f is a single component.");
    return scale_expansion(e_begin, e_end, f, h_begin, h_end);
}

template
<
    std::size_t e_length,
    std::size_t f_length,
    bool inplace,
    typename InIter,
    typename Real,
    typename OutIter,
    std::size_t result = 2 * f_length
>
inline OutIter expansion_times(Real e,
                               InIter f_begin,
                               InIter f_end,
                               OutIter h_begin,
                               OutIter h_end)
{
    static_assert(e_length == 1, "e_length must be 1 if e is a single component.");
    return scale_expansion(f_begin, f_end, e, h_begin, h_end);
}

template
<
    std::size_t e_length,
    std::size_t f_length,
    bool inplace,
    typename Real,
    typename OutIter,
    std::size_t result = 2
>
inline OutIter expansion_times(Real e,
                               Real f,
                               OutIter h_begin,
                               OutIter h_end)
{
    static_assert(e_length == 1 && f_length == 1, "e_length and f_length must be 1 if they are single components.");
    *(h_begin + 1) = e * f;
    *(h_begin) = two_product_tail(e, f, *(h_begin + 1));
    return h_begin + 2;
}

template <std::size_t s1, std::size_t s2>
struct expansion_sum_length
{
    static constexpr std::size_t value = s1 + s2;
};

template<std::size_t s1, std::size_t s2>
struct expansion_product_length
{
    static constexpr std::size_t value = 2 * s1 * s2;
};

template
<
    typename Expression,
    bool StageB = false,
    operator_types Op = Expression::operator_type
>
struct expansion_size_impl {};

template <typename Expression, bool StageB>
struct expansion_size_impl<Expression, StageB, operator_types::no_op>
{
    static constexpr std::size_t value = 1;
};

template <typename Expression, bool StageB>
struct expansion_size_impl<Expression, StageB, operator_types::sum>
{
private:
    static constexpr std::size_t left_size =
        expansion_size_impl<typename Expression::left, StageB>::value;
    static constexpr std::size_t right_size =
        expansion_size_impl<typename Expression::right, StageB>::value;
public:
    static constexpr std::size_t value = left_size + right_size;
};

template <typename Expression, bool StageB>
struct expansion_size_impl<Expression, StageB, operator_types::difference>
{
private:
    static constexpr std::size_t left_size =
        expansion_size_impl<typename Expression::left, StageB>::value;
    static constexpr std::size_t right_size =
        expansion_size_impl<typename Expression::right, StageB>::value;
public:
    static constexpr std::size_t value =
        StageB && Expression::left::is_leaf && Expression::right::is_leaf ?
              1
            : left_size + right_size;
};

template <typename Expression, bool StageB>
struct expansion_size_impl<Expression, StageB, operator_types::product>
{
private:
    static constexpr std::size_t left_size =
        expansion_size_impl<typename Expression::left, StageB>::value;
    static constexpr std::size_t right_size =
        expansion_size_impl<typename Expression::right, StageB>::value;
public:
    static constexpr std::size_t value = 2 * left_size * right_size;
};

template <typename Expression> using expansion_size =
    boost::mp11::mp_size_t< expansion_size_impl<Expression, false>::value >;

template <typename Expression> using expansion_size_stage_b =
    boost::mp11::mp_size_t< expansion_size_impl<Expression, true>::value >;

template
<
    operator_types Op,
    std::size_t LeftLength,
    std::size_t RightLength,
    bool Inplace,
    typename Iter,
    bool StageB = false
>
struct perform_op_impl {};

template
<
    std::size_t LeftLength,
    std::size_t RightLength,
    bool Inplace,
    typename Iter,
    bool StageB
>
struct perform_op_impl
    <
        operator_types::sum,
        LeftLength,
        RightLength,
        Inplace,
        Iter,
        StageB
    >
{
    template<typename ...Args>
    static Iter apply(Args...args)
    {
        return expansion_plus<LeftLength, RightLength, Inplace>(args...);
    }
};

template
<
    std::size_t LeftLength,
    std::size_t RightLength,
    bool Inplace,
    typename Iter,
    bool StageB
>
struct perform_op_impl
    <
        operator_types::difference,
        LeftLength,
        RightLength,
        Inplace,
        Iter,
        StageB
    >
{
    template<typename ...Args>
    static Iter apply(Args...args)
    {
        return expansion_minus
            <
                LeftLength,
                RightLength,
                Inplace,
                StageB
            >(args...);
    }
};

template
<
    std::size_t LeftLength,
    std::size_t RightLength,
    bool Inplace,
    typename Iter,
    bool StageB
>
struct perform_op_impl
    <
        operator_types::product,
        LeftLength,
        RightLength,
        Inplace,
        Iter,
        StageB
    >
{
    template<typename ...Args>
    static Iter apply(Args...args)
    {
        return expansion_times<LeftLength, RightLength, Inplace>(args...);
    }
};

template
<
    typename Evals,
    typename Eval,
    typename Sizes,
    typename AccumulatedSizes,
    typename Iter,
    typename Real,
    bool StageB = false,
    operator_types Op = Eval::operator_type,
    bool LeftLeaf = Eval::left::is_leaf,
    bool RightLeaf = Eval::right::is_leaf
>
struct eval_expansion_impl {};

template
<
    typename Evals,
    typename Eval,
    typename Sizes,
    typename AccumulatedSizes,
    typename Iter,
    typename Real,
    bool StageB,
    operator_types Op
>
struct eval_expansion_impl
    <
        Evals,
        Eval,
        Sizes,
        AccumulatedSizes,
        Iter,
        Real,
        StageB,
        Op,
        true,
        true
    >
{
private:
    using left = typename Eval::left;
    using right = typename Eval::right;
    using eval_index = boost::mp11::mp_find<Evals, Eval>;
    static constexpr std::size_t size =
        boost::mp11::mp_at<Sizes, eval_index>::value;
    static constexpr std::size_t start =
        boost::mp11::mp_at<AccumulatedSizes, eval_index>::value;
public:
    template<typename ...Reals>
    static Iter apply(Iter begin, Iter end, const Reals&... args)
    {
        std::array<Real, sizeof...(Reals)> input
            {{ static_cast<Real>(args)... }};
        Real left_val = input[left::argn];
        Real right_val = input[right::argn];
        return perform_op_impl<Op, 1, 1, false, Iter, StageB>
            ::apply(left_val, right_val, begin + start, begin + start + size);
    }
};

template
<
    typename Evals,
    typename Eval,
    typename Sizes,
    typename AccumulatedSizes,
    typename Iter,
    typename Real,
    bool StageB,
    operator_types Op
>
struct eval_expansion_impl
    <
        Evals,
        Eval,
        Sizes,
        AccumulatedSizes,
        Iter,
        Real,
        StageB,
        Op,
        true,
        false
    >
{
private:
    using left = typename Eval::left;
    using right = typename Eval::right;
    using eval_index = boost::mp11::mp_find<Evals, Eval>;
    static constexpr std::size_t size =
        boost::mp11::mp_at<Sizes, eval_index>::value;
    static constexpr std::size_t start =
        boost::mp11::mp_at<AccumulatedSizes, eval_index>::value;
    using right_eval_index = boost::mp11::mp_find<Evals, right>;
    static constexpr std::size_t right_size =
        boost::mp11::mp_at<Sizes, right_eval_index>::value;
    static constexpr std::size_t right_start =
        boost::mp11::mp_at<AccumulatedSizes, right_eval_index>::value;
public:
    template<typename ...Reals>
    static Iter apply(Iter begin, Iter end, const Reals&... args)
    {
        std::array<Real, sizeof...(Reals)> input
            {{ static_cast<Real>(args)... }};
        Real left_val = input[left::argn];
        return perform_op_impl<Op, 1, right_size, false, Iter, StageB>::apply(
            left_val,
            begin + right_start,
            begin + right_start + right_size,
            begin + start,
            begin + start + size);
    }
};

template
<
    typename Evals,
    typename Eval,
    typename Sizes,
    typename AccumulatedSizes,
    typename Iter,
    typename Real,
    bool StageB,
    operator_types Op
>
struct eval_expansion_impl
    <
        Evals,
        Eval,
        Sizes,
        AccumulatedSizes,
        Iter,
        Real,
        StageB,
        Op,
        false,
        true
    >
{
private:
    using left = typename Eval::left;
    using right = typename Eval::right;
    using eval_index = boost::mp11::mp_find<Evals, Eval>;
    static constexpr std::size_t size =
        boost::mp11::mp_at<Sizes, eval_index>::value;
    static constexpr std::size_t start =
        boost::mp11::mp_at<AccumulatedSizes, eval_index>::value;
    using left_eval_index = boost::mp11::mp_find<Evals, left>;
    static constexpr std::size_t left_size =
        boost::mp11::mp_at<Sizes, left_eval_index>::value;
    static constexpr std::size_t left_start =
        boost::mp11::mp_at<AccumulatedSizes, left_eval_index>::value;
public:
    template<typename ...Reals>
    static Iter apply(Iter begin, Iter end, const Reals&... args)
    {
        std::array<Real, sizeof...(Reals)> input
            {{ static_cast<Real>(args)... }};
        Real right_val = input[right::argn];
        return perform_op_impl<Op, left_size, 1, false, Iter, StageB>::apply(
            begin + left_start,
            begin + left_start + left_size,
            right_val,
            begin + start,
            begin + start + size);
    }
};

template
<
    typename Evals,
    typename Eval,
    typename Sizes,
    typename AccumulatedSizes,
    typename Iter,
    typename Real,
    bool StageB,
    operator_types Op
>
struct eval_expansion_impl
    <
        Evals,
        Eval,
        Sizes,
        AccumulatedSizes,
        Iter,
        Real,
        StageB,
        Op,
        false,
        false
    >
{
private:
    using left = typename Eval::left;
    using right = typename Eval::right;
    using eval_index = boost::mp11::mp_find<Evals, Eval>;
    static constexpr std::size_t size =
        boost::mp11::mp_at<Sizes, eval_index>::value;
    static constexpr std::size_t start =
        boost::mp11::mp_at<AccumulatedSizes, eval_index>::value;
    using left_eval_index =
        boost::mp11::mp_find<Evals, left>;
    static constexpr std::size_t left_size =
        boost::mp11::mp_at<Sizes, left_eval_index>::value;
    static constexpr std::size_t left_start =
        boost::mp11::mp_at<AccumulatedSizes, left_eval_index>::value;
    using right_eval_index = boost::mp11::mp_find<Evals, right>;
    static constexpr std::size_t right_size =
        boost::mp11::mp_at<Sizes, right_eval_index>::value;
    static constexpr std::size_t right_start =
        boost::mp11::mp_at<AccumulatedSizes, right_eval_index>::value;
public:
    template<typename ...Reals>
    static Iter apply(Iter begin, Iter end, const Reals&... args)
    {
        return perform_op_impl
            <
                Op,
                left_size,
                right_size,
                false,
                Iter,
                StageB
            >::apply(
                begin + left_start,
                begin + left_start + left_size,
                begin + right_start,
                begin + right_start + right_size,
                begin + start,
                begin + start + size);
    }
};

template
<
    typename Evals,
    typename RemainingEvals,
    typename Sizes,
    typename AccumulatedSizes,
    typename Iter,
    typename Real,
    bool StageB = false,
    std::size_t RemainingSize = boost::mp11::mp_size<RemainingEvals>::value
>
struct eval_expansions_impl
{
    template<typename ...Reals>
    static Iter apply(Iter begin, Iter end, const Reals&... args)
    {
        eval_expansion_impl
            <
                Evals,
                boost::mp11::mp_front<RemainingEvals>,
                Sizes,
                AccumulatedSizes,
                Iter,
                Real,
                StageB
            >::apply(begin, end, args...);
        return eval_expansions_impl
            <
                Evals,
                boost::mp11::mp_pop_front<RemainingEvals>,
                Sizes,
                AccumulatedSizes,
                Iter,
                Real,
                StageB
            >::apply(begin, end, args...);
    }
};

template
<
    typename Evals,
    typename RemainingEvals,
    typename Sizes,
    typename AccumulatedSizes,
    typename Iter,
    typename Real,
    bool StageB
>
struct eval_expansions_impl
    <
        Evals,
        RemainingEvals,
        Sizes,
        AccumulatedSizes,
        Iter,
        Real,
        StageB,
        1
    >
{
    template<typename ...Reals>
    static Iter apply(Iter begin, Iter end, const Reals&... args)
    {
        return eval_expansion_impl
            <
                Evals,
                boost::mp11::mp_front<RemainingEvals>,
                Sizes,
                AccumulatedSizes,
                Iter,
                Real,
                StageB
            >::apply(begin, end, args...);
    }
};

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPANSION_ARITHMETIC_HPP
