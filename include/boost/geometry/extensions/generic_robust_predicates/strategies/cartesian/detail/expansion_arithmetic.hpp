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
#include <algorithm>
#include <iterator>
#include <type_traits>

#include <boost/mp11/integral.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/function.hpp>
#include <boost/mp11/utility.hpp>
#include <boost/mp11/algorithm.hpp>

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
                std::frexp(lesser, &exp_lesser);
                std::frexp(*it, &exp_now);
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
    bool zero_elimination,
    typename OutIter,
    typename Real
>
struct insert_ze_impl
{
    static constexpr OutIter apply(OutIter out, Real val)
    {
        if(val == Real(0))
        {
            return out;
        }
        else
        {
            *out = val;
            ++out;
            return out;
        }
    }
};

template
<
    typename OutIter,
    typename Real
>
struct insert_ze_impl<false, OutIter, Real>
{
    static constexpr OutIter apply(OutIter out, Real val)
    {
        *out = val;
        ++out;
        return out;
    }
};

template
<
    bool zero_elimination,
    typename OutIter,
    typename Real
>
constexpr OutIter insert_ze(OutIter o, Real r)
{
    return insert_ze_impl<zero_elimination, OutIter, Real>::apply(o, r);
}

template
<
    bool zero_elimination,
    typename OutIter,
    typename Real
>
struct insert_ze_final_impl
{
    static constexpr OutIter apply(OutIter out, OutIter start, Real val)
    {
        if(val == Real(0) && out != start)
        {
            return out;
        }
        else
        {
            *out = val;
            ++out;
            return out;
        }
    }
};

template
<
    typename OutIter,
    typename Real
>
struct insert_ze_final_impl<false, OutIter, Real>
{
    static constexpr OutIter apply(OutIter out, OutIter, Real val)
    {
        *out = val;
        ++out;
        return out;
    }
};

template
<
    bool zero_elimination,
    typename OutIter,
    typename Real
>
constexpr OutIter insert_ze_final(OutIter o, OutIter s, Real r)
{
    return insert_ze_final_impl<zero_elimination, OutIter, Real>::apply(o, s, r);
}

template
<
    bool zero_elimination,
    typename InIter,
    typename OutIter,
    typename Real,
    bool NegateE = false,
    bool NegateB = false
>
constexpr OutIter grow_expansion(InIter e_begin,
                                 InIter e_end,
                                 Real b,
                                 OutIter h_begin,
                                 OutIter)
{
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));
    Real Q = negate<NegateB>(b);
    auto h_it = h_begin;
    for(auto e_it = e_begin; e_it != e_end; ++e_it)
    {
        Real Q_new = negate<NegateE>(*e_it) + Q;
        Real h_new = two_sum_tail(negate<NegateE>(*e_it), Q, Q_new);
        Q = Q_new;
        h_it = insert_ze<zero_elimination>(h_it, h_new);
    }
    h_it = insert_ze_final<zero_elimination>(h_it, h_begin, Q);
    assert(debug_expansion::expansion_nonoverlapping(h_begin, h_it));
    assert(  !debug_expansion::expansion_nonadjacent(e_begin, e_end)
           || debug_expansion::expansion_nonadjacent(h_begin, h_it) );
    return h_it;
}

template
<
    bool zero_elimination,
    typename InIter,
    typename OutIter,
    typename Real
>
constexpr OutIter grow_expansion_difference(InIter e_begin,
                                            InIter e_end,
                                            Real b,
                                            OutIter h_begin,
                                            OutIter h_end)
{
    return grow_expansion
        <
            zero_elimination,
            InIter,
            OutIter,
            Real,
            false,
            true
        >(e_begin, e_end, b, h_begin, h_end);
}

template <bool zero_elimination, bool NegateE, bool NegateB>
struct expansion_sum_advance_impl
{
    template <typename Real>
    static constexpr bool apply(Real, Real)
    {
        return true;
    }
};

template <bool NegateE, bool NegateB>
struct expansion_sum_advance_impl<true, NegateE, NegateB>
{
    template <typename Real>
    static constexpr bool apply(Real e, Real b)
    {
        Real Q = negate<NegateE>(e) + negate<NegateB>(b);
        return two_sum_tail(negate<NegateE>(e), negate<NegateB>(b), Q) != Real(0);
    }
};


template <bool zero_elimination, bool NegateE, bool NegateB, typename Real>
constexpr bool expansion_sum_advance(Real e, Real b)
{
    return expansion_sum_advance_impl<zero_elimination, NegateE, NegateB>
            ::apply(e, b);
}

template
<
    bool zero_elimination,
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool NegateE = false,
    bool NegateF = false
>
constexpr OutIter expansion_sum(InIter1 e_begin,
                                InIter1 e_end,
                                InIter2 f_begin,
                                InIter2 f_end,
                                OutIter h_begin,
                                OutIter)
{
    using Real = typename std::iterator_traits<InIter1>::value_type;
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));
    assert(debug_expansion::expansion_nonoverlapping(f_begin, f_end));
    const auto elen = std::distance(e_begin, e_end);
    const auto flen = std::distance(f_begin, f_end);
    auto f_it = f_begin;
    auto h_begin_i = h_begin;
    auto h_it = grow_expansion
        <
            zero_elimination,
            InIter1,
            OutIter,
            Real,
            NegateE,
            NegateF
        >(e_begin, e_end, *f_it, h_begin, h_begin + elen + 1);
    if(expansion_sum_advance<zero_elimination, NegateE, NegateF>(*e_begin,
                                                                 *f_it))
    {
        ++h_begin_i;
    }
    ++f_it;
    for(auto i = 1; i < flen; ++i)
    {
        h_it = grow_expansion
            <
                zero_elimination,
                InIter1,
                OutIter,
                Real,
                false,
                NegateF
            >(h_begin_i,
              h_it,
              *f_it,
              h_begin_i,
              h_it + 1);

        if(expansion_sum_advance<zero_elimination, false, NegateF>(*h_begin_i,
                                                                   *f_it))
        {
            ++h_begin_i;
        }
        ++f_it;
    }
    assert(debug_expansion::expansion_nonoverlapping(h_begin, h_it));
    assert(   !debug_expansion::expansion_nonadjacent(e_begin, e_end)
           || !debug_expansion::expansion_nonadjacent(f_begin, f_end)
           || debug_expansion::expansion_nonadjacent(h_begin, h_it));
    return h_it;
}

template
<
    bool zero_elimination,
    typename InIter1,
    typename InIter2,
    typename OutIter
>
constexpr OutIter expansion_difference(InIter1 e_begin,
                                       InIter1 e_end,
                                       InIter2 f_begin,
                                       InIter2 f_end,
                                       OutIter h_begin,
                                       OutIter h_end)
{
    return expansion_sum
        <
            zero_elimination,
            InIter1,
            InIter2,
            OutIter,
            false,
            true
        >(e_begin, e_end, f_begin, f_end, h_begin, h_end);
}

template
<
    bool zero_elimination,
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool NegateE = false,
    bool NegateF = false,
    bool e_no_zeros = false,
    bool f_no_zeros = false
>
constexpr OutIter fast_expansion_sum_inplace(InIter1 e_begin,
                                             InIter1 e_end,
                                             InIter2 f_begin,
                                             InIter2 f_end,
                                             OutIter h_begin,
                                             OutIter h_end)
{
    assert(e_end == f_begin);
    assert(f_begin != h_begin);
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

    auto e_old_end = e_end;
    typename std::iterator_traits<InIter1>::difference_type e_shortened = 0;
    if(!e_no_zeros)
    {
        e_end = std::remove(e_begin, e_end, Real(0));
        e_shortened = std::distance(e_end, e_old_end);
    }
    std::rotate(e_end, f_begin, f_end);
    f_begin = e_end;
    auto f_old_end = f_end;
    f_end = f_end - e_shortened;
    if(!f_no_zeros)
    {
        f_end = std::remove(f_begin, f_end, Real(0));
    }
    if(!zero_elimination)
    {
        std::fill(f_end, f_old_end, Real(0));
    }

    std::inplace_merge(e_begin, e_end, f_end, abs_comp{});

    InIter1 g_it = e_begin;
    InIter2 g_end = f_end;
    auto h_it = h_begin;
    Real Q = *g_it + *(g_it + 1);
    auto h_new = fast_two_sum_tail(*(g_it + 1), *(g_it), Q);
    h_it = insert_ze<zero_elimination>(h_it, h_new);
    g_it += 2;
    for(; g_it < g_end; ++g_it)
    {
        Real Q_new = Q + *g_it;
        h_new = two_sum_tail(Q, *g_it, Q_new);
        h_it = insert_ze<zero_elimination>(h_it, h_new);
        Q = Q_new;
    }
    h_it = insert_ze_final<zero_elimination>(h_it, h_begin, Q);
    //assert(debug_expansion::expansion_strongly_nonoverlapping(h_begin, h_it));
    if(zero_elimination)
    {
        return h_it;
    }
    else
    {
        return h_end;
    }
}

template
<
    bool zero_elimination,
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool NegateE = false,
    bool NegateF = false
>
constexpr OutIter fast_expansion_sum_not_inplace(InIter1 e_begin,
                                                 InIter1 e_end,
                                                 InIter2 f_begin,
                                                 InIter2 f_end,
                                                 OutIter h_begin,
                                                 OutIter)
{
    assert(e_begin != h_begin);
    assert(f_begin != h_begin);
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
        h_it = insert_ze<zero_elimination>(h_it, h_new);
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
            h_it = insert_ze<zero_elimination>(h_it, h_new);
        }
    }
    while(e_it != e_end)
    {
        Real Q_new = negate<NegateE>(*e_it) + Q;
        Real h_new = two_sum_tail(negate<NegateE>(*e_it), Q, Q_new);
        h_it = insert_ze<zero_elimination>(h_it, h_new);
        Q = Q_new;
        ++e_it;
    }
    while(f_it != f_end)
    {
        Real Q_new = negate<NegateF>(*f_it) + Q;
        Real h_new = two_sum_tail(negate<NegateF>(*f_it), Q, Q_new);
        h_it = insert_ze<zero_elimination>(h_it, h_new);
        Q = Q_new;
        ++f_it;
    }
    h_it = insert_ze_final<zero_elimination>(h_it, h_begin, Q);
    //assert(debug_expansion::expansion_strongly_nonoverlapping(h_begin, h_it));
    return h_it;
}

template
<
    bool zero_elimination,
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool inplace,
    bool NegateE,
    bool NegateF
>
struct fast_expansion_sum_impl
{
    static constexpr OutIter apply(InIter1 e_begin,
                                   InIter1 e_end,
                                   InIter2 f_begin,
                                   InIter2 f_end,
                                   OutIter h_begin,
                                   OutIter h_end)
    {
        return fast_expansion_sum_not_inplace
            <
                zero_elimination,
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
    bool zero_elimination,
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool NegateE,
    bool NegateF
>
struct fast_expansion_sum_impl
    <
        zero_elimination,
        InIter1,
        InIter2,
        OutIter,
        true,
        NegateE,
        NegateF
    >
{
    static constexpr OutIter apply(InIter1 e_begin,
                                   InIter1 e_end,
                                   InIter2 f_begin,
                                   InIter2 f_end,
                                   OutIter h_begin,
                                   OutIter h_end)
    {
        return fast_expansion_sum_inplace
            <
                zero_elimination,
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
    bool zero_elimination,
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool inplace,
    bool NegateE,
    bool NegateF
>
constexpr OutIter fast_expansion_sum(InIter1 e_begin,
                                     InIter1 e_end,
                                     InIter2 f_begin,
                                     InIter2 f_end,
                                     OutIter h_begin,
                                     OutIter h_end)
{
    return fast_expansion_sum_impl
        <
            zero_elimination,
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
    bool zero_elimination,
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
            zero_elimination,
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
    bool zero_elimination,
    typename InIter,
    typename Real,
    typename OutIter
>
constexpr OutIter scale_expansion(InIter e_begin,
                                  InIter e_end,
                                  Real b,
                                  OutIter h_begin,
                                  OutIter)
{
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));

    auto e_it = e_begin;
    auto h_it = h_begin;
    Real Q = *e_it * b;
    auto h_new = two_product_tail(*e_it, b, Q);
    h_it = insert_ze<zero_elimination>(h_it, h_new);
    ++e_it;
    for(; e_it != e_end; ++e_it)
    {
        Real product_1 = *e_it * b;
        Real product_0 = two_product_tail(*e_it, b, product_1);
        Real sum = Q + product_0;
        h_new = two_sum_tail(Q, product_0, sum);
        h_it = insert_ze<zero_elimination>(h_it, h_new);
        Q = product_1 + sum;
        h_new = two_sum_tail(product_1, sum, Q);
        h_it = insert_ze<zero_elimination>(h_it, h_new);
    }
    h_it = insert_ze_final<zero_elimination>(h_it, h_begin, Q);

    assert( debug_expansion::expansion_nonoverlapping(h_begin, h_it) );
    assert(  !debug_expansion::expansion_nonadjacent(e_begin, e_end)
           || debug_expansion::expansion_nonadjacent(h_begin, h_it) );
    assert(  !debug_expansion::expansion_strongly_nonoverlapping(e_begin, e_end)
           || debug_expansion::expansion_strongly_nonoverlapping(h_begin, h_it) );
    return h_it;
}

template <typename Rhs>
struct greater_than_or_equal
{
    template<typename Lhs>
    using fn = boost::mp11::mp_bool< Lhs::value >= Rhs::value >;
};

template <int s1, int s2>
struct expansion_sum_length
{
    static constexpr int value = s1 != -1 && s2 != -1 ? s1 + s2 : -1;
};

template<int s1, int s2>
struct expansion_product_length
{
    static constexpr int value = s1 != -1 && s2 != -1 ? 2 * s1 * s2 : -1;
};

template
<
    int e_length,
    int f_length,
    bool inplace = false,
    bool e_negate = false,
    bool f_negate = false,
    int h_length = expansion_sum_length<e_length, f_length>::value
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
                (false),
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
struct expansion_plus_impl<1, 1, inplace, e_negate, f_negate>
{
    template<typename InIter1, typename InIter2, typename OutIter>
    static inline OutIter apply(InIter1 e_begin,
                                InIter1,
                                InIter2 f_begin,
                                InIter2,
                                OutIter h_begin,
                                OutIter)
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
    int e_length,
    bool inplace,
    bool e_negate,
    bool f_negate,
    int h_length
>
struct expansion_plus_impl<e_length, 1, inplace, e_negate, f_negate, h_length>
{
    template<typename InIter1, typename InIter2, typename OutIter>
    static inline OutIter apply(InIter1 e_begin,
                                InIter1 e_end,
                                InIter2 f_begin,
                                InIter2,
                                OutIter h_begin,
                                OutIter h_end)
    {
        return grow_expansion
            <
                (false),
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
    int f_length,
    bool inplace,
    bool e_negate,
    bool f_negate,
    int h_length
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
                                InIter1,
                                InIter2 f_begin,
                                InIter2 f_end,
                                OutIter h_begin,
                                OutIter h_end)
    {
        return grow_expansion
            <
                (false),
                InIter2,
                OutIter,
                typename std::iterator_traits<InIter1>::value_type,
                e_negate,
                f_negate
            >(f_begin, f_end, *e_begin, h_begin, h_end);
    }
};

template<bool inplace, bool e_negate, bool f_negate>
struct expansion_plus_impl<2, 2, inplace, e_negate, f_negate>
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
                false,
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
    int e_length,
    int f_length,
    bool inplace,
    typename InIter1,
    typename InIter2,
    typename OutIter,
    int result = expansion_sum_length<e_length, f_length>::value
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
    int e_length,
    int f_length,
    bool inplace,
    typename InIter,
    typename Real,
    typename OutIter,
    int result = expansion_sum_length<e_length, f_length>::value
>
inline OutIter expansion_plus(InIter e_begin,
                              InIter e_end,
                              Real f,
                              OutIter h_begin,
                              OutIter h_end)
{
    static_assert( f_length == 1, "f_length must be 1 if f is a single component." );
    return grow_expansion<(false)>(e_begin, e_end, f, h_begin, h_end);
}

template
<
    int e_length,
    int f_length,
    bool inplace,
    typename InIter,
    typename Real,
    typename OutIter,
    int result = expansion_sum_length<e_length, f_length>::value
>
inline OutIter expansion_plus(Real e,
                              InIter f_begin,
                              InIter f_end,
                              OutIter h_begin,
                              OutIter h_end)
{
    static_assert( e_length == 1, "e_length must be 1 if e is a single component." );
    return grow_expansion<(false)>(f_begin, f_end, e, h_begin, h_end);
}

template
<
    int e_length,
    int f_length,
    bool inplace,
    typename Real,
    typename OutIter,
    int result = expansion_sum_length<e_length, f_length>::value
>
inline OutIter expansion_plus(
    Real e,
    Real f,
    OutIter h_begin,
    OutIter)
{
    static_assert( f_length == 1 && e_length == 1, "e_length and f_length must be 1 if they are single components." );
    *(h_begin + 1) = e + f;
    *h_begin = two_sum_tail(e, f, *(h_begin + 1));
    return h_begin + 2;
}

template
<
    int e_length,
    int f_length,
    bool inplace,
    bool StageB,
    typename InIter1,
    typename InIter2,
    typename OutIter,
    int result = expansion_sum_length<e_length, f_length>::value
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
    int e_length,
    int f_length,
    bool inplace,
    bool StageB,
    typename InIter,
    typename Real,
    typename OutIter,
    int result = expansion_sum_length<e_length, f_length>::value
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
    int e_length,
    int f_length,
    bool inplace,
    bool StageB,
    typename Real,
    typename InIter,
    typename OutIter,
    int result = expansion_sum_length<e_length, f_length>::value
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
            (false),
            InIter,
            OutIter,
            Real,
            false,
            true
        >(f_begin, f_end, e, h_begin, h_end);
}

template
<
    int e_length,
    int f_length,
    bool inplace,
    bool StageB,
    typename Real,
    typename OutIter,
    int result = expansion_sum_length<e_length, f_length>::value
>
inline std::enable_if_t<!StageB, OutIter> expansion_minus(Real e,
                                                         Real f,
                                                         OutIter h_begin,
                                                         OutIter)
{
    static_assert(e_length == 1 && f_length == 1, "e_length and f_length must be 1 if they are single components.");
    *(h_begin + 1) = e - f;
    *h_begin = two_difference_tail(e, f, *(h_begin + 1));
    return h_begin + 2;
}

template
<
    int e_length,
    int f_length,
    bool inplace,
    bool StageB,
    typename Real,
    typename OutIter,
    int result = expansion_sum_length<e_length, f_length>::value
>
inline std::enable_if_t<StageB, OutIter> expansion_minus(Real e,
                                                         Real f,
                                                         OutIter h_begin,
                                                         OutIter)
{
    static_assert(e_length == 1 && f_length == 1, "e_length and f_length must be 1 if they are single components.");
    *(h_begin) = e - f;
    return h_begin + 1;
}

template
<
    int e_length,
    int f_length,
    bool inplace,
    typename InIter,
    typename Real,
    typename OutIter,
    int result = expansion_product_length<e_length, f_length>::value
>
inline OutIter expansion_times(InIter e_begin,
                               InIter e_end,
                               Real f,
                               OutIter h_begin,
                               OutIter h_end)
{
    static_assert(f_length == 1, "f_length must be 1 if f is a single component.");
    return scale_expansion<(false)>(e_begin, e_end, f, h_begin, h_end);
}

template
<
    int e_length,
    int f_length,
    bool inplace,
    typename InIter,
    typename Real,
    typename OutIter,
    int result = expansion_product_length<e_length, f_length>::value
>
inline OutIter expansion_times(Real e,
                               InIter f_begin,
                               InIter f_end,
                               OutIter h_begin,
                               OutIter h_end)
{
    static_assert(e_length == 1, "e_length must be 1 if e is a single component.");
    return scale_expansion<(false)>(f_begin, f_end, e, h_begin, h_end);
}

template
<
    int e_length,
    int f_length,
    bool inplace,
    typename Real,
    typename OutIter,
    int result = expansion_product_length<e_length, f_length>::value
>
inline OutIter expansion_times(Real e,
                               Real f,
                               OutIter h_begin,
                               OutIter)
{
    static_assert(e_length == 1 && f_length == 1, "e_length and f_length must be 1 if they are single components.");
    *(h_begin + 1) = e * f;
    *(h_begin) = two_product_tail(e, f, *(h_begin + 1));
    return h_begin + 2;
}

template
<
    int e_length,
    int f_length,
    bool,
    typename In1,
    typename In2,
    typename Out,
    int result = expansion_product_length<e_length, f_length>::value
>
constexpr Out expansion_times(In1, In1, In2, In2, Out, Out);

template
<
    int e_length,
    int f_length,
    int result = expansion_product_length<e_length, f_length>::value,
    bool e_smaller = e_length <= f_length
>
struct expansion_times_impl
{
    template<typename InIter1, typename InIter2, typename OutIter>
    static constexpr OutIter apply(InIter1 e_begin,
                                   InIter1 e_end,
                                   InIter2 f_begin,
                                   InIter2 f_end,
                                   OutIter h_begin,
                                   OutIter h_end)
    {
        assert(e_begin != h_begin && f_begin != h_begin);

        //TODO: Evaluate zero-elimination for very short expansions before multiplication.
        const auto e_dyn_length = std::distance(e_begin, e_end);
        if(e_dyn_length == 1)
        {
            return expansion_times<1, f_length, false>(*e_begin,
                                                       f_begin,
                                                       f_end,
                                                       h_begin,
                                                       h_end);
        }
        else if (e_dyn_length > 1 )
        {
            constexpr int e_length1 = e_length == -1 ? -1 : e_length / 2;
            auto e_mid = e_begin + e_dyn_length / 2;
            auto h_mid = expansion_times
                <
                    e_length1,
                    f_length,
                    false
                >(e_begin, e_mid, f_begin, f_end, h_begin, h_end);

            constexpr int e_length2 = e_length == -1 ? -1 : e_length - e_length1;
            h_end = expansion_times
                <
                    e_length2,
                    f_length,
                    false
                >(e_mid, e_end, f_begin, f_end, h_mid, h_end);

            constexpr int summand_length1 =
                expansion_product_length<e_length1, f_length>::value;
            constexpr int summand_length2 =
                expansion_product_length<e_length2, f_length>::value;
            auto h_it = expansion_plus
                <
                    summand_length1,
                    summand_length2,
                    true
                >(h_begin, h_mid, h_mid, h_end, h_begin, h_end);
            assert(debug_expansion::expansion_nonoverlapping(h_begin, h_it));
            return h_it;
        }
        else if (e_dyn_length == 0)
        {
            return h_begin;
        }
        assert(false);
    }
};

template
<
    int e_length,
    int f_length,
    int result
>
struct expansion_times_impl<e_length, f_length, result, false>
{
    template<typename InIter1, typename InIter2, typename OutIter>
    static constexpr OutIter apply(InIter1 e_begin,
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
    static constexpr OutIter apply(InIter1 e_begin,
                                   InIter1,
                                   InIter2 f_begin,
                                   InIter2,
                                   OutIter h_begin,
                                   OutIter)
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
    int e_length,
    int f_length,
    bool inplace,
    typename InIter1,
    typename InIter2,
    typename OutIter,
    int result = expansion_product_length<e_length, f_length>::value
>
constexpr OutIter expansion_times(InIter1 e_begin,
                                  InIter1 e_end,
                                  InIter2 f_begin,
                                  InIter2 f_end,
                                  OutIter h_begin,
                                  OutIter h_end)
{
    return expansion_times_impl<e_length, f_length>::
        apply(e_begin, e_end, f_begin, f_end, h_begin, h_end);
}

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPANSION_ARITHMETIC_HPP
