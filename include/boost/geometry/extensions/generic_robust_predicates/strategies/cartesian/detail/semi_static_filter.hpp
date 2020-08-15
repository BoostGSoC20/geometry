// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SEMI_STATIC_FILTER_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SEMI_STATIC_FILTER_HPP

#include <array>

#include <boost/mp11/list.hpp>
#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/set.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/approximate.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/result_propagation.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template
<
    typename Expression,
    typename CalculationType,
    typename ErrorExpression
>
struct semi_static_filter
{
private:
    using stack = typename boost::mp11::mp_unique<post_order<Expression>>;
    using evals = typename boost::mp11::mp_remove_if<stack, is_leaf>;
    using error_eval_stack = boost::mp11::mp_unique
        <
            post_order<ErrorExpression>
        >;
    using error_eval_stack_remainder = boost::mp11::mp_set_difference
        <
            error_eval_stack,
            evals
        >;
    using all_evals = boost::mp11::mp_append
        <
            evals,
            error_eval_stack_remainder
        >;
    using ct = CalculationType;
    static constexpr std::size_t expr_max_arg = max_argn<Expression>::value;
    static constexpr std::size_t error_expr_max_arg =
        max_argn<ErrorExpression>::value;
public:
    using computations = boost::mp11::mp_list<Expression, ErrorExpression>;
    static constexpr bool stateful = false;
    static constexpr bool updates = false;
    static constexpr std::size_t arg_count =
        expr_max_arg > error_expr_max_arg ? expr_max_arg : error_expr_max_arg;

    template
    <
        typename InputList,
        typename StatefulStages,
        typename RemainingStages,
        typename RemainingReusables,
        typename RemainingComputations,
        typename ...InputArr
    >
    static inline int staged_apply(const StatefulStages& stages, const InputArr&... inputs)
    {
        using reusables = boost::mp11::mp_front<RemainingReusables>;
        using comps = boost::mp11::mp_front<RemainingComputations>;
        std::array<ct, boost::mp11::mp_size<reusables>::value> reusable;
        using nonreusables = boost::mp11::mp_set_difference<comps, reusables>;
        std::array<ct, boost::mp11::mp_size<nonreusables>::value> nonreusable;
        using arg_list = boost::mp11::mp_push_back
            <
                InputList,
                reusables,
                nonreusables
            >;
        approximate_interim<comps, arg_list, ct>(inputs...,
                                                 reusable,
                                                 nonreusable);
        const ct error_bound =
            get_approx<ErrorExpression, arg_list, ct>(inputs...,
                                                      reusable,
                                                      nonreusable);
        const ct det = get_approx<Expression, arg_list, ct>(inputs...,
                                                            reusable,
                                                            nonreusable);
        if (det > error_bound)
        {
            return 1;
        }
        else if (det < -error_bound)
        {
            return -1;
        }
        else if (error_bound == 0 && det == 0)
        {
            return 0;
        }
        else
        {
            using next_input_list =
                boost::mp11::mp_push_back<InputList, reusables>;
            using next_remaining_stages =
                boost::mp11::mp_pop_front<RemainingStages>;
            using next_remaining_reusables =
                boost::mp11::mp_pop_front<RemainingReusables>;
            using next_remaining_computations =
                boost::mp11::mp_pop_front<RemainingComputations>;
            return next_stage
                <
                    next_input_list,
                    StatefulStages,
                    next_remaining_stages,
                    next_remaining_reusables,
                    next_remaining_computations
                >::apply(stages, inputs..., reusable);
        }
    }

    template <typename ...Reals>
    static inline int apply(const Reals&... args)
    {
        using arg_list_input = argument_list<sizeof...(Reals)>;
        using arg_list = boost::mp11::mp_list<all_evals, arg_list_input>;
        std::array<CalculationType, sizeof...(Reals)> input
            {{ static_cast<ct>(args)... }};
        std::array<ct, boost::mp11::mp_size<all_evals>::value> results;
        approximate_interim<all_evals, arg_list, ct>(results, input);
        const ct error_bound =
            get_approx<ErrorExpression, arg_list, ct>(results, input);
        const ct det = get_approx<Expression, arg_list, ct>(results, input);
        if (det > error_bound)
        {
            return 1;
        }
        else if (det < -error_bound)
        {
            return -1;
        }
        else if (error_bound == 0 && det == 0)
        {
            return 0;
        }
        else
        {
            return sign_uncertain;
        }
    }
};

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SEMI_STATIC_FILTER_HPP
