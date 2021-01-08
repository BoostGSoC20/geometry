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

constexpr unsigned long long round_to_power_of_two(unsigned long long num)
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
    return y;
}

template <typename Real>
constexpr Real fast_two_sum_tail(Real a, Real b, Real x)
{
    assert(std::abs(a) >= std::abs(b) || a == 0);
    Real b_virtual = x - a;
    Real y = b - b_virtual;
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
    return y;
}

template <typename Real>
constexpr Real fast_two_difference_tail(Real a, Real b, Real x)
{
    assert(std::abs(a) >= std::abs(b) || a == 0);
    Real b_virtual = a - x;
    Real y = b_virtual - b;
    return y;
}

template <typename Real>
constexpr Real two_product_tail(Real a, Real b, Real x)
{
    Real y = std::fma(a, b, -x);
    return y;
}

template <typename Real>
constexpr Real splitter_helper()
{
    int digits = std::numeric_limits<Real>::digits;
    int half = digits / 2;
    int ceilhalf = half + (half * 2 == digits ? 0 : 1 );
    Real out(1);
    for(int i = 0; i < ceilhalf; ++i)
    {
        out *= Real(2);
    }
    return out + Real(1);
}

template <typename Real>
constexpr Real splitter = splitter_helper<Real>();

template <typename Real>
constexpr std::array<Real, 2> split(Real a)
{
    Real c = splitter<Real> * a;
    Real a_big = c - a;
    Real a_hi = c - a_big;
    Real a_lo = a - a_hi;
    return std::array<Real, 2>{ a_hi, a_lo };
}

template <typename Real>
constexpr Real two_product_tail_constexpr(Real a, Real b, Real x)
{
    const auto a_split = split(a);
    const auto b_split = split(b);
    Real a_hi_b_hi = a_split[0] * b_split[0];
    Real err1 = x - a_hi_b_hi;
    Real a_lo_b_hi = a_split[1] * b_split[0];
    Real err2 = err1 - a_lo_b_hi;
    Real a_hi_b_lo = a_split[0] * b_split[1];
    Real err3 = err2 - a_hi_b_lo;
    Real a_lo_b_lo = a_split[1] * b_split[1];
    return a_lo_b_lo - err3;
}

template
<
    bool ZeroElimination,
    bool MostSigOnly,
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
struct insert_ze_impl<false, true, OutIter, Real>
{
    static constexpr OutIter apply(OutIter out, Real val)
    {
        *out += val;
        return out;
    }
};

template
<
    typename OutIter,
    typename Real
>
struct insert_ze_impl<true, true, OutIter, Real>
{
    static constexpr OutIter apply(OutIter out, Real val)
    {
        if(val != 0)
        {
            *out = val;
        }
        return out;
    }
};

template
<
    typename OutIter,
    typename Real
>
struct insert_ze_impl<false, false, OutIter, Real>
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
    bool ZeroElimination,
    bool MostSigOnly = false,
    typename OutIter,
    typename Real
>
constexpr OutIter insert_ze(OutIter o, Real r)
{
    return insert_ze_impl<ZeroElimination, MostSigOnly, OutIter, Real>::apply(o, r);
}

template
<
    bool ZeroElimination,
    bool MostSigOnly,
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
struct insert_ze_final_impl<false, true, OutIter, Real>
{
    static constexpr OutIter apply(OutIter out, OutIter, Real val)
    {
        *out += val;
        ++out;
        return out;
    }
};

template
<
    typename OutIter,
    typename Real
>
struct insert_ze_final_impl<true, true, OutIter, Real>
{
    static constexpr OutIter apply(OutIter out, OutIter, Real val)
    {
        if(val != 0)
        {
            *out = val;
        }
        ++out;
        return out;
    }
};

template
<
    typename OutIter,
    typename Real
>
struct insert_ze_final_impl<false, false, OutIter, Real>
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
    bool ZeroElimination,
    bool MostSigOnly = false,
    typename OutIter,
    typename Real
>
constexpr OutIter insert_ze_final(OutIter o, OutIter s, Real r)
{
    return insert_ze_final_impl<ZeroElimination, MostSigOnly, OutIter, Real>::apply(o, s, r);
}

template
<
    bool ZeroElimination,
    bool MostSigOnly = false,
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
        h_it = insert_ze<ZeroElimination>(h_it, h_new);
    }
    h_it = insert_ze_final<ZeroElimination>(h_it, h_begin, Q);
    assert(debug_expansion::expansion_nonoverlapping(h_begin, h_it));
    assert(  !debug_expansion::expansion_nonadjacent(e_begin, e_end)
           || debug_expansion::expansion_nonadjacent(h_begin, h_it) );
    return h_it;
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
    bool ZeroElimination,
    bool MostSigOnly = false,
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
    bool advance =
        expansion_sum_advance<ZeroElimination, NegateE, NegateF>(*e_begin,
                                                                 *f_it);
    auto h_it = grow_expansion
        <
            ZeroElimination,
            false,
            InIter1,
            OutIter,
            Real,
            NegateE,
            NegateF
        >(e_begin, e_end, *f_it, h_begin, h_begin + elen + 1);
    if(advance)
    {
        ++h_begin_i;
    }
    ++f_it;
    for(auto i = 1; i < flen; ++i)
    {
        advance =
            expansion_sum_advance<ZeroElimination, false, NegateF>(*h_begin_i,
                                                                    *f_it);
        h_it = grow_expansion
            <
                ZeroElimination,
                false,
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
        if(advance)
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
    bool ZeroElimination,
    bool MostSigOnly = false,
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
    assert(e_end <= f_begin);
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
    auto f_old_end = f_end;
    f_end = f_end - e_shortened - std::distance(e_old_end, f_begin);
    f_begin = e_end;

    if(!f_no_zeros)
    {
        f_end = std::remove(f_begin, f_end, Real(0));
    }
    if(!ZeroElimination)
    {
        std::fill(f_end, f_old_end, Real(0));
    }

    std::inplace_merge(e_begin, e_end, f_end, abs_comp{});

    InIter1 g_it = e_begin;
    InIter2 g_end = f_end;
    OutIter h_it = h_begin;
    Real Q = *g_it + *(g_it + 1);
    Real h_new = fast_two_sum_tail(*(g_it + 1), *(g_it), Q);
    h_it = insert_ze<ZeroElimination, MostSigOnly>(h_it, h_new);
    g_it += 2;
    for(; g_it < g_end; ++g_it)
    {
        Real Q_new = Q + *g_it;
        h_new = two_sum_tail(Q, *g_it, Q_new);
        h_it = insert_ze<ZeroElimination, MostSigOnly>(h_it, h_new);
        Q = Q_new;
    }
    h_it = insert_ze_final<ZeroElimination, MostSigOnly>(h_it, h_begin, Q);

    assert(debug_expansion::expansion_nonoverlapping(h_begin, h_it));
    if( !ZeroElimination )
    {
        h_it = h_end;
    }
    return h_it;
}

template
<
    bool ZeroElimination,
    bool MostSigOnly = false,
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
                                                 OutIter h_end)
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
        h_it = insert_ze<ZeroElimination, MostSigOnly>(h_it, h_new);
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
            h_it = insert_ze<ZeroElimination, MostSigOnly>(h_it, h_new);
        }
    }
    while(e_it != e_end)
    {
        Real Q_new = negate<NegateE>(*e_it) + Q;
        Real h_new = two_sum_tail(negate<NegateE>(*e_it), Q, Q_new);
        h_it = insert_ze<ZeroElimination, MostSigOnly>(h_it, h_new);
        Q = Q_new;
        ++e_it;
    }
    while(f_it != f_end)
    {
        Real Q_new = negate<NegateF>(*f_it) + Q;
        Real h_new = two_sum_tail(negate<NegateF>(*f_it), Q, Q_new);
        h_it = insert_ze<ZeroElimination, MostSigOnly>(h_it, h_new);
        Q = Q_new;
        ++f_it;
    }
    h_it = insert_ze_final<ZeroElimination, MostSigOnly>(h_it, h_begin, Q);
    assert(debug_expansion::expansion_nonoverlapping(h_begin, h_it));
    return h_it;
}

template
<
    bool ZeroElimination,
    bool MostSigOnly,
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
                ZeroElimination,
                MostSigOnly,
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
    bool ZeroElimination,
    bool MostSigOnly,
    typename InIter1,
    typename InIter2,
    typename OutIter,
    bool NegateE,
    bool NegateF
>
struct fast_expansion_sum_impl
    <
        ZeroElimination,
        MostSigOnly,
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
                ZeroElimination,
                MostSigOnly,
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
    bool ZeroElimination,
    bool MostSigOnly = false,
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
            ZeroElimination,
            MostSigOnly,
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

constexpr int expansion_sum_length(int s1, int s2)
{
    return s1 != -1 && s2 != -1 ? s1 + s2 : -1;
}

constexpr int expansion_product_length(int s1, int s2, bool same = false)
{
    if(s1 == -1 || s2 == -1)
    {
        return -1;
    }
    else if(same && s1 == 2 && s2 == 2)
    {
        return 6;
    }
    return 2 * s1 * s2;
}

template <int Length>
using default_zero_elimination_policy =
    boost::mp11::mp_bool<(Length > 16) || Length == -1>;

template <int Length1, int Length2>
using default_fast_expansion_sum_policy =
    boost::mp11::mp_bool<   (Length1 > 2 || Length1 == -1)
                         && (Length2 > 2 || Length2 == -2)>;

template
<
    int e_length,
    int f_length,
    bool inplace = false,
    bool e_negate = false,
    bool f_negate = false,
    template<int> class zero_elimination = default_zero_elimination_policy,
    bool fast_expansion = default_fast_expansion_sum_policy
        <
            e_length,
            f_length
        >::value,
    bool MostSigOnly = false,
    int h_length = expansion_sum_length(e_length, f_length)
>
struct expansion_plus_impl
{
    template<typename InIter1, typename InIter2, typename OutIter>
    static constexpr OutIter apply(
        InIter1 e_begin,
        InIter1 e_end,
        InIter2 f_begin,
        InIter2 f_end,
        OutIter h_begin,
        OutIter h_end)
    {
        return fast_expansion_sum
            <
                zero_elimination<h_length>::value,
                MostSigOnly,
                InIter1,
                InIter2,
                OutIter,
                inplace,
                e_negate,
                f_negate
            >(e_begin, e_end, f_begin, f_end, h_begin, h_end);
    }
};

template <bool inplace, bool e_negate, bool f_negate, template<int> class ze, bool MostSigOnly, int h_length>
struct expansion_plus_impl<1, 1, inplace, e_negate, f_negate, ze, false, MostSigOnly, h_length>
{
    template<typename InIter1, typename InIter2, typename OutIter>
    static constexpr OutIter apply(InIter1 e_begin,
                                   InIter1,
                                   InIter2 f_begin,
                                   InIter2,
                                   OutIter h_begin,
                                   OutIter)
    {
        auto x = negate<e_negate>(*e_begin) + negate<f_negate>(*f_begin);
        auto y = two_sum_tail(negate<e_negate>(*e_begin),
                              negate<f_negate>(*f_begin), x);
        auto h_it = h_begin;
        h_it = insert_ze<ze<h_length>::value>(h_it, y);
        h_it = insert_ze_final<ze<h_length>::value>(h_it, h_begin, x);
        return h_it;
    }
};

template
<
    int e_length,
    bool inplace,
    bool e_negate,
    bool f_negate,
    template<int> class ze,
    bool MostSigOnly,
    int h_length
>
struct expansion_plus_impl<e_length, 1, inplace, e_negate, f_negate, ze, false, MostSigOnly, h_length>
{
    template<typename InIter1, typename InIter2, typename OutIter>
    static constexpr OutIter apply(InIter1 e_begin,
                                   InIter1 e_end,
                                   InIter2 f_begin,
                                   InIter2,
                                   OutIter h_begin,
                                   OutIter h_end)
    {
        return grow_expansion
            <
                ze<h_length>::value,
                MostSigOnly,
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
    template<int> class ze,
    bool fe,
    bool MostSigOnly,
    int h_length
>
struct expansion_plus_impl
    <
        1,
        f_length,
        inplace,
        e_negate,
        f_negate,
        ze,
        fe,
        MostSigOnly,
        h_length
    >
{
    template<typename InIter1, typename InIter2, typename OutIter>
    static constexpr OutIter apply(InIter1 e_begin,
                                   InIter1,
                                   InIter2 f_begin,
                                   InIter2 f_end,
                                   OutIter h_begin,
                                   OutIter h_end)
    {
        return grow_expansion
            <
                ze<h_length>::value,
                MostSigOnly,
                InIter2,
                OutIter,
                typename std::iterator_traits<InIter1>::value_type,
                e_negate,
                f_negate
            >(f_begin, f_end, *e_begin, h_begin, h_end);
    }
};

template
<
    int f_length,
    int e_length,
    bool inplace,
    bool e_negate,
    bool f_negate,
    template<int> class ze,
    bool MostSigOnly,
    int h_length
>
struct expansion_plus_impl<f_length, e_length, inplace, e_negate, f_negate, ze, false, MostSigOnly, h_length>
{
    template<typename InIter1, typename InIter2, typename OutIter>
    static constexpr OutIter apply(InIter1 e_begin,
                                   InIter1 e_end,
                                   InIter2 f_begin,
                                   InIter2 f_end,
                                   OutIter h_begin,
                                   OutIter h_end)
    {
        return expansion_sum
            <
                ze<h_length>::value,
                MostSigOnly,
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
    template<int> class ze = default_zero_elimination_policy,
    template<int, int> class fe = default_fast_expansion_sum_policy,
    bool MostSigOnly = false,
    typename InIter1,
    typename InIter2,
    typename OutIter,
    int result = expansion_sum_length(e_length, f_length)
>
constexpr OutIter expansion_plus(InIter1 e_begin,
                                 InIter1 e_end,
                                 InIter2 f_begin,
                                 InIter2 f_end,
                                 OutIter h_begin,
                                 OutIter h_end)
{
    return expansion_plus_impl
        <
            e_length,
            f_length,
            inplace,
            false,
            false,
            ze,
            fe<e_length, f_length>::value,
            MostSigOnly
        >::apply(e_begin, e_end, f_begin, f_end, h_begin, h_end);
}

template
<
    int e_length,
    int f_length,
    bool inplace,
    template<int> class ze = default_zero_elimination_policy,
    template<int, int> class = default_fast_expansion_sum_policy,
    bool MostSigOnly = false,
    typename InIter,
    typename Real,
    typename OutIter,
    int result = expansion_sum_length(e_length, f_length)
>
constexpr OutIter expansion_plus(InIter e_begin,
                                 InIter e_end,
                                 Real f,
                                 OutIter h_begin,
                                 OutIter h_end)
{
    static_assert( f_length == 1, "f_length must be 1 if f is a single component." );
    return grow_expansion<ze<result>::value, MostSigOnly>(e_begin, e_end, f, h_begin, h_end);
}

template
<
    int e_length,
    int f_length,
    bool inplace,
    template<int> class ze = default_zero_elimination_policy,
    template<int, int> class = default_fast_expansion_sum_policy,
    bool MostSigOnly = false,
    typename InIter,
    typename Real,
    typename OutIter,
    int result = expansion_sum_length(e_length, f_length)
>
constexpr OutIter expansion_plus(Real e,
                                 InIter f_begin,
                                 InIter f_end,
                                 OutIter h_begin,
                                 OutIter h_end)
{
    static_assert( e_length == 1, "e_length must be 1 if e is a single component." );
    return grow_expansion<ze<result>::value, MostSigOnly>(f_begin, f_end, e, h_begin, h_end);
}

template
<
    int e_length,
    int f_length,
    bool inplace,
    template<int> class ze = default_zero_elimination_policy,
    template<int, int> class = default_fast_expansion_sum_policy,
    bool MostSigOnly = false,
    typename Real,
    typename OutIter,
    int result = expansion_sum_length(e_length, f_length)
>
constexpr OutIter expansion_plus(
    Real e,
    Real f,
    OutIter h_begin,
    OutIter)
{
    static_assert( f_length == 1 && e_length == 1, "e_length and f_length must be 1 if they are single components." );
    Real x = e + f;
    if(MostSigOnly)
    {
        *h_begin = x;
        return h_begin + 1;
    }
    Real y = two_sum_tail(e, f, *(h_begin + 1));
    OutIter h_it = h_begin;
    h_it = insert_ze<ze<result>::value>(h_it, y);
    h_it = insert_ze_final<ze<result>::value>(h_it, x);
    return h_it;
}

template
<
    int e_length,
    int f_length,
    bool inplace,
    bool StageB,
    template<int> class ze = default_zero_elimination_policy,
    template<int, int> class fe = default_fast_expansion_sum_policy,
    bool MostSigOnly = false,
    typename InIter1,
    typename InIter2,
    typename OutIter,
    int result = expansion_sum_length(e_length, f_length)
>
constexpr OutIter expansion_minus(InIter1 e_begin,
                                  InIter1 e_end,
                                  InIter2 f_begin,
                                  InIter2 f_end,
                                  OutIter h_begin,
                                  OutIter h_end)
{
    return expansion_plus_impl
        <
            e_length,
            f_length,
            inplace,
            false,
            true,
            ze,
            fe<e_length, f_length>::value,
            MostSigOnly
        >::apply(e_begin, e_end, f_begin, f_end, h_begin, h_end);
}

template
<
    int e_length,
    int f_length,
    bool inplace,
    bool StageB,
    template<int> class ZE = default_zero_elimination_policy,
    template<int, int> class FE = default_fast_expansion_sum_policy,
    bool MostSigOnly = false,
    typename InIter,
    typename Real,
    typename OutIter,
    int result = expansion_sum_length(e_length, f_length)
>
constexpr OutIter expansion_minus(InIter e_begin,
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
            inplace,
            false,
            true,
            ZE,
            FE,
            MostSigOnly
        >(e_begin, e_end, f, h_begin, h_end);
}

template
<
    int e_length,
    int f_length,
    bool inplace,
    bool StageB,
    template<int> class ze = default_zero_elimination_policy,
    template<int, int> class = default_fast_expansion_sum_policy,
    bool MostSigOnly = false,
    typename Real,
    typename InIter,
    typename OutIter,
    int result = expansion_sum_length(e_length, f_length)
>
constexpr OutIter expansion_minus(Real e,
                                  InIter f_begin,
                                  InIter f_end,
                                  OutIter h_begin,
                                  OutIter h_end)
{
    static_assert(e_length == 1, "e_length must be 1 if e is a single component.");
    return grow_expansion
        <
            ze<result>::value,
            MostSigOnly,
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
    template<int> class ze = default_zero_elimination_policy,
    template<int, int> class = default_fast_expansion_sum_policy,
    bool MostSigOnly = false,
    typename Real,
    typename OutIter,
    int result = expansion_sum_length(e_length, f_length)
>
constexpr std::enable_if_t<!StageB, OutIter> expansion_minus(Real e,
                                                             Real f,
                                                             OutIter h_begin,
                                                             OutIter)
{
    static_assert(e_length == 1 && f_length == 1, "e_length and f_length must be 1 if they are single components.");
    Real x = e - f;
    if (MostSigOnly)
    {
        *h_begin = x;
        return h_begin + 1;
    }
    Real y = two_difference_tail(e, f, x);
    OutIter h_it = h_begin;
    h_it = insert_ze<ze<result>::value>(h_it, y);
    h_it = insert_ze_final<ze<result>::value>(h_it, h_begin, x);
    return h_it;
}

template
<
    int e_length,
    int f_length,
    bool inplace,
    bool StageB,
    template<int> class = default_zero_elimination_policy,
    template<int, int> class = default_fast_expansion_sum_policy,
    bool MostSigOnly = false,
    typename Real,
    typename OutIter,
    int result = expansion_sum_length(e_length, f_length)
>
constexpr std::enable_if_t<StageB, OutIter> expansion_minus(Real e,
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
    template<int> class ze = default_zero_elimination_policy,
    template<int, int> class = default_fast_expansion_sum_policy,
    bool LeftEqualsRight,
    typename InIter,
    typename Real,
    typename OutIter,
    int result = expansion_product_length(e_length, f_length, LeftEqualsRight)
>
constexpr OutIter expansion_times(InIter e_begin,
                                  InIter e_end,
                                  Real f,
                                  OutIter h_begin,
                                  OutIter h_end)
{
    static_assert(f_length == 1, "f_length must be 1 if f is a single component.");
    return scale_expansion<ze<result>::value>(e_begin, e_end, f, h_begin, h_end);
}

template
<
    int e_length,
    int f_length,
    bool inplace,
    template<int> class ze = default_zero_elimination_policy,
    template<int, int> class = default_fast_expansion_sum_policy,
    bool LeftEqualsRight,
    typename InIter,
    typename Real,
    typename OutIter,
    int result = expansion_product_length(e_length, f_length, LeftEqualsRight)
>
constexpr OutIter expansion_times(Real e,
                                  InIter f_begin,
                                  InIter f_end,
                                  OutIter h_begin,
                                  OutIter h_end)
{
    static_assert(e_length == 1, "e_length must be 1 if e is a single component.");
    auto h_it = scale_expansion<ze<result>::value>(f_begin, f_end, e, h_begin, h_end);
    return h_it;
}

template
<
    int e_length,
    int f_length,
    bool inplace,
    template<int> class ze = default_zero_elimination_policy,
    template<int, int> class = default_fast_expansion_sum_policy,
    bool LeftEqualsRight,
    typename Real,
    typename OutIter,
    int result = expansion_product_length(e_length, f_length, LeftEqualsRight)
>
constexpr OutIter expansion_times(Real e,
                                  Real f,
                                  OutIter h_begin,
                                  OutIter)
{
    static_assert(e_length == 1 && f_length == 1, "e_length and f_length must be 1 if they are single components.");
    Real x = e * f;
    Real y = two_product_tail(e, f, x);
    OutIter h_it = h_begin;
    h_it = insert_ze<ze<result>::value>(h_it, y);
    h_it = insert_ze_final<ze<result>::value>(h_it, h_begin, x);
    return h_it;
}

template
<
    int e_length,
    int f_length,
    bool,
    template<int> class = default_zero_elimination_policy,
    template<int, int> class = default_fast_expansion_sum_policy,
    bool LeftEqualsRight,
    typename In1,
    typename In2,
    typename Out,
    int = expansion_product_length(e_length, f_length, LeftEqualsRight)
>
constexpr Out expansion_times(In1, In1, In2, In2, Out, Out);

template
<
    int e_length,
    int f_length,
    template <int> class ZeroEliminationPolicy =
        default_zero_elimination_policy,
    template<int, int> class FastExpansionSumPolicy =
        default_fast_expansion_sum_policy,
    bool LeftEqualsRight = false,
    int result = expansion_product_length(e_length, f_length, LeftEqualsRight),
    typename e_smaller = boost::mp11::mp_bool<e_length <= f_length>
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
        assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));
        assert(debug_expansion::expansion_nonoverlapping(f_begin, f_end));
        //TODO: Evaluate zero-elimination for very short expansions before multiplication.
        const auto e_dyn_length = std::distance(e_begin, e_end);
        const auto h_old_end = h_end;
        if(e_dyn_length == 1)
        {
            auto h_it = expansion_times
                <
                    1,
                    f_length,
                    false,
                    ZeroEliminationPolicy,
                    FastExpansionSumPolicy,
                    LeftEqualsRight
                >(*e_begin, f_begin, f_end, h_begin, h_end);
            assert(debug_expansion::expansion_nonoverlapping(h_begin, h_it));
            return h_it;
        }
        else if (e_dyn_length > 1 )
        {
            constexpr int e_length1 = e_length == -1 ? -1 : e_length / 2;
            auto e_mid = e_begin + e_dyn_length / 2;
            auto h_mid = h_begin + std::distance(e_begin, e_mid) * std::distance(f_begin, f_end) * 2;
            h_mid = expansion_times
                <
                    e_length1,
                    f_length,
                    false,
                    ZeroEliminationPolicy,
                    FastExpansionSumPolicy,
                    LeftEqualsRight
                >(e_begin, e_mid, f_begin, f_end, h_begin, h_mid);
            constexpr int e_length2 = e_length == -1 ? -1 : e_length - e_length1;
            h_end = expansion_times
                <
                    e_length2,
                    f_length,
                    false,
                    ZeroEliminationPolicy,
                    FastExpansionSumPolicy,
                    LeftEqualsRight
                >(e_mid, e_end, f_begin, f_end, h_mid, h_end);

            constexpr int summand_length1 =
                expansion_product_length(e_length1, f_length, false);
            constexpr int summand_length2 =
                expansion_product_length(e_length2, f_length, false);
            OutIter h_it = expansion_plus
                <
                    summand_length1,
                    summand_length2,
                    true,
                    ZeroEliminationPolicy,
                    FastExpansionSumPolicy
                >(h_begin, h_mid, h_mid, h_end, h_begin, h_end);
            assert(debug_expansion::expansion_nonoverlapping(h_begin, h_it));
            return h_it;
        }
        else if (e_dyn_length == 0)
        {
            return h_begin;
        }
        assert(false);
        return h_begin;
    }
};

/*
template
<
    int e_length,
    int f_length,
    bool inplace,
    template<int> class ze = default_zero_elimination_policy,
    template<int, int> class = default_fast_expansion_sum_policy,
    typename InIter,
    typename Real,
    typename OutIter,
    int result = expansion_sum_length(e_length, f_length)
>
constexpr OutIter expansion_plus(InIter e_begin,
                                 InIter e_end,
                                 Real f,
                                 OutIter h_begin,
                                 OutIter h_end)

 * */

template
<
    bool ZeroElimination
>
struct advance_ze_impl
{
    template <typename Iter>
    Iter apply(Iter it)
    {
        return *it == 0 ? it : it + 1;
    }
};

template <>
struct advance_ze_impl<false>
{
    template <typename Iter>
    Iter apply(Iter it)
    {
        return it + 1;
    }
};

template <bool B>
struct constant_ze
{
    template <int>
    using fn = boost::mp11::mp_bool<B>;
};

template
<
    bool ze,
    bool EZeroEliminated
>
struct two_square_impl
{
    template <typename InIter, typename OutIter>
    static constexpr OutIter apply(InIter e_begin,
                                   InIter e_end,
                                   OutIter h_begin,
                                   OutIter h_end)
    {
        std::array<std::remove_reference_t<decltype(*h_begin)>, 5> cache;
        cache[2] = *(e_begin) * *(e_begin);
        auto h_it = h_begin;
        h_it = insert_ze<ze>(h_it,
                             two_product_tail(*e_begin, *e_begin, cache[2]));
        cache[1] = *(e_begin + 1) * 2.0 * *(e_begin);
        cache[0] = two_product_tail(*(e_begin + 1), 2.0 * *(e_begin), cache[1]);
        expansion_plus<2, 1, false>(cache.cbegin(),
                                    cache.cbegin() + 2,
                                    cache[2],
                                    cache.begin() + 2,
                                    cache.begin() + 5);
        h_it = insert_ze<ze>(h_it, cache[2]);
        cache[1] = *(e_begin + 1) * *(e_begin + 1);
        cache[0] = two_product_tail(*(e_begin + 1), *(e_begin + 1), cache[1]);
        h_it = expansion_plus
            <
                2,
                2,
                false,
                constant_ze<ze>::template fn
            >(cache.cbegin(),
              cache.cbegin() + 2,
              cache.cbegin() + 3,
              cache.cbegin() + 5,
              h_it,
              h_it + 4);
        return h_it;
    }
};

template <>
struct two_square_impl<true, true>
{
    template <typename InIter, typename OutIter>
    static constexpr OutIter apply(InIter e_begin,
                                   InIter e_end,
                                   OutIter h_begin,
                                   OutIter h_end)
    {
        if(e_begin + 1 == e_end)
        {
            *(h_begin + 1) = *e_begin * *e_begin;
            *h_begin = two_product_tail(*e_begin, *e_begin, *(h_begin + 1));
            return h_begin + 2;
        }
        else
        {
            return two_square_impl<true, false>
                ::apply(e_begin, e_end, h_begin, h_end);
        }
    }
};

template
<
    template<int> class ZE,
    template<int, int> class FE
>
struct expansion_times_impl
    <
        2,
        2,
        ZE,
        FE,
        true,
        6,
        boost::mp11::mp_true
    >
{
    template<typename InIter1, typename InIter2, typename OutIter>
    static constexpr OutIter apply(InIter1 e_begin,
                                   InIter1 e_end,
                                   InIter2,
                                   InIter2,
                                   OutIter h_begin,
                                   OutIter h_end)
    {
        return two_square_impl
            <
                ZE<6>::value,
                ZE<2>::value
            >::apply(e_begin, e_end, h_begin, h_end);
    }
};


template
<
    int e_length,
    int f_length,
    template<int> class ze,
    template<int, int> class fe,
    bool LeftEqualsRight,
    int result
>
struct expansion_times_impl
    <
        e_length,
        f_length,
        ze,
        fe,
        LeftEqualsRight,
        result,
        boost::mp11::mp_false
    >
{
    template<typename InIter1, typename InIter2, typename OutIter>
    static constexpr OutIter apply(InIter1 e_begin,
                                   InIter1 e_end,
                                   InIter2 f_begin,
                                   InIter2 f_end,
                                   OutIter h_begin,
                                   OutIter h_end)
    {
        return expansion_times_impl<f_length, e_length, ze, fe, LeftEqualsRight>
            ::apply(f_begin,
                    f_end,
                    e_begin,
                    e_end,
                    h_begin,
                    h_end);
    }
};

template<template<int> class ze, template<int, int> class fe, bool LeftEqualsRight>
struct expansion_times_impl<1, 1, ze, fe, LeftEqualsRight, 2>
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
        OutIter h_it = h_begin;
        h_it = insert_ze<ze<2>::value>(h_it, y);
        h_it = insert_ze_final<ze<2>::value>(h_it, h_begin, x);
        return h_it;
    }
};

template
<
    int e_length,
    int f_length,
    bool inplace,
    template<int> class ze,
    template<int, int> class fe,
    bool LeftEqualsRight,
    typename InIter1,
    typename InIter2,
    typename OutIter,
    int result
>
constexpr OutIter expansion_times(InIter1 e_begin,
                                  InIter1 e_end,
                                  InIter2 f_begin,
                                  InIter2 f_end,
                                  OutIter h_begin,
                                  OutIter h_end)
{
    return expansion_times_impl<e_length, f_length, ze, fe, LeftEqualsRight, result>::
        apply(e_begin, e_end, f_begin, f_end, h_begin, h_end);
}

template
<
    typename Iter
>
constexpr Iter compress(Iter e_begin, Iter e_end)
{
    Iter bottom_it = std::prev(e_end);
    auto Q = *bottom_it;
    Iter e_it;
    for (e_it = std::prev(bottom_it); e_it != e_begin; --e_it)
    {
        auto Q_next = Q + *e_it;
        auto q = fast_two_sum_tail(Q, *e_it, Q_next);
        Q = Q_next;
        if (q != 0)
        {
            *bottom_it = Q;
            --bottom_it;
            Q = q;
        }
    }
    auto Q_next = Q + *e_it;
    auto q = fast_two_sum_tail(Q, *e_it, Q_next);
    Q = Q_next;
    if (q != 0)
    {
        *bottom_it = Q;
        --bottom_it;
        Q = q;
    }
    *bottom_it = Q;
    Iter top_it = e_begin;
    for(e_it = std::next(bottom_it); e_it != e_end; ++e_it)
    {
        auto Q_next = *e_it + Q;
        auto q = fast_two_sum_tail(*e_it, Q, Q_next);
        Q = Q_next;
        if(q != 0)
        {
            *top_it = q;
            ++top_it;
        }
    }
    *top_it = Q;
    ++top_it;
    return top_it;
}

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPANSION_ARITHMETIC_HPP