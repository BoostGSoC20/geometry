// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGE_B_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGE_B_HPP

#include <array>
#include <iterator>

#include <boost/mp11/list.hpp>
#include <boost/mp11/algorithm.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_d.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template
<
    typename Expression,
    bool IsDifference = Expression::operator_type == operator_types::difference
>
struct is_leaf_difference_impl
{
    using type = boost::mp11::mp_false;
};

template
<
    typename Expression
>
struct is_leaf_difference_impl <Expression, true>
{
    using type = boost::mp11::mp_bool
        <
            Expression::left::is_leaf && Expression::right::is_leaf
        >;
};

template <typename Expression>
using is_leaf_difference = typename is_leaf_difference_impl<Expression>::type;

template
<
    typename Evals,
    typename LeafDifferences,
    typename AccumulatedSizes,
    typename Iter,
    typename Real,
    std::size_t RemainingDifferences =
        boost::mp11::mp_size<LeafDifferences>::value
>
struct all_differences_zero_tail
{
    template <typename ...Reals>
    static bool apply(Iter begin, Iter end, const Reals&... args)
    {
        using eval = boost::mp11::mp_front<LeafDifferences>;
        using left = typename eval::left;
        using right = typename eval::right;
        std::array<Real, sizeof...(Reals)> input
            {{ static_cast<Real>(args)... }};
        Real left_val = input[left::argn - 1];
        Real right_val = input[right::argn - 1];
        using eval_index = boost::mp11::mp_find<Evals, eval>;
        constexpr std::size_t start =
            boost::mp11::mp_at<AccumulatedSizes, eval_index>::value;
        if(two_difference_tail(left_val, right_val, *(begin + start))
                == Real(0))
        {
            return all_differences_zero_tail
                <
                    Evals,
                    boost::mp11::mp_pop_front<LeafDifferences>,
                    AccumulatedSizes,
                    Iter,
                    Real
                >::apply(begin, end, args...);
        }
        return false;
    }
};

template
<
    typename Evals,
    typename LeafDifferences,
    typename AccumulatedSizes,
    typename Iter,
    typename Real
>
struct all_differences_zero_tail
    <
        Evals,
        LeafDifferences,
        AccumulatedSizes,
        Iter,
        Real,
        0
    >
{
    template <typename ...Reals>
    static bool apply(Iter, Iter, const Reals&...)
    {
        return true;
    }
};

template <typename Expression, typename Real>
struct stage_b
{
    static constexpr bool stateful = false;
    static constexpr bool updates = false;

    template <typename ...Reals>
    static inline int apply(const Reals&... args)
    {
        using stack = typename boost::mp11::mp_unique<post_order<Expression>>;
        using evals = typename boost::mp11::mp_remove_if<stack, is_leaf>;
        using sizes = boost::mp11::mp_transform<expansion_size_stage_b, evals>;
        using accumulated_sizes = boost::mp11::mp_push_front
            <
                boost::mp11::mp_partial_sum
                    <
                        sizes,
                        boost::mp11::mp_size_t<0>,
                        boost::mp11::mp_plus
                    >,
                boost::mp11::mp_size_t<0>
            >;

        using result_array =
            std::array<Real, boost::mp11::mp_back<accumulated_sizes>::value>;
        result_array results;

        using leaf_differences =
            boost::mp11::mp_copy_if<evals, is_leaf_difference>;

        /*auto final_exp_end = */eval_expansions_impl
            <
                evals,
                evals,
                sizes,
                accumulated_sizes,
                decltype(results.begin()),
                Real,
                true
            >::apply(results.begin(), results.end(), args...);

        bool all_zero = all_differences_zero_tail
                <
                    evals,
                    leaf_differences,
                    accumulated_sizes,
                    typename result_array::iterator,
                    Real
                >::apply(results.begin(), results.end(), args...);
        if( !all_zero )
        {
            return sign_uncertain;
        }

        constexpr std::size_t final_exp_size =
            boost::mp11::mp_back<sizes>::value;
        auto is_zero = [](Real d) { return d == Real(0); };
        auto most_significant = std::find_if_not(
                results.crbegin(),
                results.crbegin() + final_exp_size, is_zero);
        if( most_significant == results.crbegin() + final_exp_size)
        {
            return 0;
        }
        else if( *most_significant > 0 )
        {
            return 1;
        }
        else
        {
            return -1;
        }
    }
};

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGE_B_HPP
