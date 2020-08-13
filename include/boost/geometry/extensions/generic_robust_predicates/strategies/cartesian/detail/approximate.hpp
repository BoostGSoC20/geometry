// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_APPROXIMATE_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_APPROXIMATE_HPP

#include <cstddef>
#include <cmath>
#include <array>
#include <tuple>

#include <boost/mp11/integral.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/algorithm.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/result_propagation.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template
<
    typename Node,
    typename InputList,
    typename Real,
    bool IsStaticConstant,
    typename ...InputArr
>
struct get_approx_impl
{
    static inline Real apply(const InputArr&... inputs)
    {
        using indices = index_pair<Node, InputList>;
        return std::get<indices::second_type::value>(
                std::get<indices::first_type::value>(
                    std::forward_as_tuple(inputs...)
                )
            );
    }
};

template
<
    typename Node,
    typename InputList,
    typename Real,
    typename ...InputArr
>
struct get_approx_impl
    <
        Node,
        InputList,
        Real,
        true,
        InputArr...
    >
{
    static inline Real apply(const InputArr&...)
    {
        return Node::value;
    }
};

template
<
    typename Node,
    typename InputList,
    typename Real,
    typename ...InputArr
>
inline Real get_approx(const InputArr&... inputs)
{
    return get_approx_impl
        <
            Node,
            InputList,
            Real,
            is_static_constant<Node>::value,
            InputArr...
        >::apply(inputs...);
}

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename Arr,
    operator_types Op,
    typename ...InputArr
>
struct approximate_interim_impl {};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename Arr,
    bool Empty,
    typename ...InputArr
>
struct approximate_remainder_impl
{
    static inline void apply(Arr& interim_results, const InputArr&... inputs)
    {
        using node = boost::mp11::mp_front<Remaining>;
        approximate_interim_impl
            <
                All,
                Remaining,
                InputList,
                Real,
                Arr,
                node::operator_type,
                InputArr...
            >::apply(interim_results, inputs...);
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename Arr,
    typename ...InputArr
>
struct approximate_remainder_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        Arr,
        true,
        InputArr...
    >
{
    static inline void apply(Arr&, const InputArr&...) {}
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename Arr,
    typename ...InputArr
>
inline void approximate_remainder(Arr& interim_results, const InputArr&... inputs)
{
    approximate_remainder_impl
        <
            All,
            Remaining,
            InputList,
            Real,
            Arr,
            boost::mp11::mp_empty<Remaining>::value,
            InputArr...
        >::apply(interim_results, inputs...);
}

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename Arr,
    typename ...InputArr
>
struct approximate_interim_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        Arr,
        operator_types::product,
        InputArr...
    >
{
    static inline void apply(Arr& interim_results, const InputArr&... inputs)
    {
        using node = boost::mp11::mp_front<Remaining>;
        using allm = boost::mp11::mp_push_front<InputList, All>;
        interim_results[boost::mp11::mp_find<All, node>::value] =
                  get_approx<typename node::left, allm, Real>(interim_results,
                                                              inputs...)
                * get_approx<typename node::right, allm, Real>(interim_results,
                                                               inputs...);
        approximate_remainder
            <
                All,
                boost::mp11::mp_pop_front<Remaining>,
                InputList,
                Real
            >(interim_results, inputs...);
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename Arr,
    typename ...InputArr
>
struct approximate_interim_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        Arr,
        operator_types::max,
        InputArr...
    >
{
    static inline void apply(Arr& interim_results, const InputArr&... inputs)
    {
        using node = boost::mp11::mp_front<Remaining>;
        using allm = boost::mp11::mp_push_front<InputList, All>;
        interim_results[boost::mp11::mp_find<All, node>::value] = std::max(
                  get_approx<typename node::left, allm, Real>(interim_results,
                                                              inputs...),
                  get_approx<typename node::right, allm, Real>(interim_results,
                                                               inputs...));
        approximate_remainder
            <
                All,
                boost::mp11::mp_pop_front<Remaining>,
                InputList,
                Real
            >(interim_results, inputs...);
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename Arr,
    typename ...InputArr
>
struct approximate_interim_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        Arr,
        operator_types::min,
        InputArr...
    >
{
    static inline void apply(Arr& interim_results, const InputArr&... inputs)
    {
        using node = boost::mp11::mp_front<Remaining>;
        using allm = boost::mp11::mp_push_front<InputList, All>;
        interim_results[boost::mp11::mp_find<All, node>::value] = std::min(
                  get_approx<typename node::left, allm, Real>(interim_results,
                                                              inputs...),
                  get_approx<typename node::right, allm, Real>(interim_results,
                                                               inputs...));
        approximate_remainder
            <
                All,
                boost::mp11::mp_pop_front<Remaining>,
                InputList,
                Real
            >(interim_results, inputs...);
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename Arr,
    typename ...InputArr
>
struct approximate_interim_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        Arr,
        operator_types::sum,
        InputArr...
    >
{
    static inline void apply(Arr& interim_results, const InputArr&... inputs)
    {
        using node = boost::mp11::mp_front<Remaining>;
        using allm = boost::mp11::mp_push_front<InputList, All>;
        interim_results[boost::mp11::mp_find<All, node>::value] =
                  get_approx<typename node::left, allm, Real>(interim_results,
                                                              inputs...)
                + get_approx<typename node::right, allm, Real>(interim_results,
                                                               inputs...);
        approximate_remainder
            <
                All,
                boost::mp11::mp_pop_front<Remaining>,
                InputList,
                Real
            >(interim_results, inputs...);
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename Arr,
    typename ...InputArr
>
struct approximate_interim_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        Arr,
        operator_types::difference,
        InputArr...
    >
{
    static inline void apply(Arr& interim_results, const InputArr&... inputs)
    {
        using node = boost::mp11::mp_front<Remaining>;
        using allm = boost::mp11::mp_push_front<InputList, All>;
        interim_results[boost::mp11::mp_find<All, node>::value] =
                  get_approx<typename node::left, allm, Real>(interim_results,
                                                              inputs...)
                - get_approx<typename node::right, allm, Real>(interim_results,
                                                               inputs...);
        approximate_remainder
            <
                All,
                boost::mp11::mp_pop_front<Remaining>,
                InputList,
                Real
            >(interim_results, inputs...);
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename Arr,
    typename ...InputArr
>
struct approximate_interim_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        Arr,
        operator_types::abs,
        InputArr...
    >
{
    static inline void apply(Arr& interim_results, const InputArr&... inputs)
    {
        using node = boost::mp11::mp_front<Remaining>;
        using allm = boost::mp11::mp_push_front<InputList, All>;
        interim_results[boost::mp11::mp_find<All, node>::value] =
            std::abs(get_approx
                <
                    typename node::child,
                    allm,
                    Real
                >(interim_results, inputs...));
        approximate_remainder
            <
                All,
                boost::mp11::mp_pop_front<Remaining>,
                InputList,
                Real
            >(interim_results, inputs...);
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename Arr,
    typename ...InputArr
>
struct approximate_interim_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        Arr,
        operator_types::no_op,
        InputArr...
    >
{
    static inline void apply(Arr& interim_results, const InputArr&... inputs)
    {
        approximate_remainder
            <
                All,
                boost::mp11::mp_pop_front<Remaining>,
                InputList,
                Real
            >(interim_results, inputs...);
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename Arr,
    typename ...InputArr
>
inline void approximate_interim(Arr& interim_results, const InputArr&... inputs)
{
    approximate_remainder
        <
            All,
            Remaining,
            InputList,
            Real
        >(interim_results, inputs...);
}

template<typename Expression, typename Real, typename InputArr>
inline Real approximate_value(const InputArr& input)
{
    using stack = typename boost::mp11::mp_unique<post_order<Expression>>;
    using evals = typename boost::mp11::mp_remove_if<stack, is_leaf>;
    using arg_list = boost::mp11::mp_list
        <
            argument_list<std::tuple_size<InputArr>::value>
        >;
    std::array<Real, boost::mp11::mp_size<evals>::value> results;
    approximate_interim<evals, evals, arg_list, Real>(results, input);
    return results.back();
}

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_APPROXIMATE_HPP
