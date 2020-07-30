// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGE_D_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGE_D_HPP

#include <array>
#include <iterator>

#include <boost/mp11/list.hpp>
#include <boost/mp11/algorithm.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/approximate.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expansion_arithmetic.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template <typename Expression, typename Real>
struct stage_d
{
    static constexpr bool stateful = false;
    static constexpr bool updates = false;
    using computations = boost::mp11::mp_list<>; //TODO: make use of previous comps

    template <typename ...Reals>
    static inline int apply(const Reals&... args)
    {
        using root = Expression;
        using stack = typename boost::mp11::mp_unique<post_order<Expression>>;
        using evals = typename boost::mp11::mp_remove_if<stack, is_leaf>;
        using sizes = boost::mp11::mp_transform<expansion_size, evals>;
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

        auto final_exp_end = eval_expansions_impl
            <
                evals,
                evals,
                sizes,
                accumulated_sizes,
                decltype(results.begin()),
                Real
            >::apply(results.begin(), results.end(), args...);


        //TODO: If the last value is not computed with zero elimination, we need to either search
        //      for the last non-zero entry (multiplate branches) or add them all up (multiple
        //      float ops). Evaluate which one is faster (if the case of no zero elimination
        //      even matters at all).
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

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGE_D_HPP
