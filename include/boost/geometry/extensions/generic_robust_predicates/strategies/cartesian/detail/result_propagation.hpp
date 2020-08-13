// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_RESULT_PROPAGATION_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_RESULT_PROPAGATION_HPP

#include <boost/mp11/integral.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/map.hpp>
#include <boost/mp11/bind.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template <typename Haystack>
struct contained_in_or_leaf_q
{
    template <typename Needle> using fn =
        boost::mp11::mp_or
            <
                boost::mp11::mp_contains<Haystack, Needle>,
                is_leaf<Needle>
            >;
};

template <typename Haystacks>
struct multi_contained_in_or_leaf_q
{
    template <typename Needle> using fn = boost::mp11::mp_or
        <
            typename multi_contained_in_or_leaf_q
                <
                    boost::mp11::mp_pop_front<Haystacks>
                >::template fn<Needle>,
            typename contained_in_or_leaf_q<boost::mp11::mp_front<Haystacks>>
                ::template fn<Needle>
        >;
};

template <>
struct multi_contained_in_or_leaf_q<boost::mp11::mp_list<>>
{
    template <typename Needle> using fn = boost::mp11::mp_false;
};

template
<
    std::size_t Last,
    typename In = boost::mp11::mp_list<>,
    std::size_t At = 1
>
struct argument_list_impl
{
private:
    using next = boost::mp11::mp_push_back<In, argument<At>>;
public:
    using type = typename argument_list_impl<Last, next, At + 1>::type;
};

template
<
    std::size_t Last,
    typename In
>
struct argument_list_impl
    <
        Last,
        In,
        Last
    >
{
    using type = boost::mp11::mp_push_back<In, argument<Last>>;
};

template <std::size_t Last> using argument_list =
    typename argument_list_impl<Last>::type;

template <template <typename> class Anchor>
struct remaining_computations_fold_q
{
    template
    <
        typename V,
        typename T
    >
    using fn = boost::mp11::mp_unique
        <
            boost::mp11::mp_append
            <
                V,
                post_order
                    <
                        T,
                        Anchor
                    >
            >
        >;
};

template
<
    typename Computations,
    typename PrevComps
>
struct remaining_computations_impl
{
    using type = boost::mp11::mp_fold_q
        <
            Computations,
            boost::mp11::mp_list<>,
            remaining_computations_fold_q
                <
                    multi_contained_in_or_leaf_q<PrevComps>::template fn
                >
        >;
};

template
<
    typename Filters,
    typename PrevComps
>
struct computations_forward {
private:
    using filter = boost::mp11::mp_front<Filters>;
    using computations = typename filter::computations;
    using remaining = typename remaining_computations_impl<computations, PrevComps>::type;
    using next_comps = boost::mp11::mp_push_back<PrevComps, remaining>;
public:
    using type = typename computations_forward<boost::mp11::mp_pop_front<Filters>, next_comps>::type;
};

template
<
    typename PrevComps
>
struct computations_forward<boost::mp11::mp_list<>, PrevComps>
{
    using type = PrevComps;
};

template
<
    typename Child,
    typename Parent
>
using is_direct_parent =
    boost::mp11::mp_contains<typename Parent::all_children, Child>;

template
<
    typename Child
>
struct parent_of_q
{
    template<typename Parent> using fn = is_direct_parent<Child, Parent>;
};

template
<
    typename Child
>
struct any_parent_of_q
{
    template<typename ParentList> using fn = boost::mp11::mp_any_of_q
        <
            ParentList,
            parent_of_q<Child>
        >;
};

template
<
    typename Expression
>
struct wanted_by_q
{
    template <typename Filter> using fn =
        boost::mp11::mp_contains<typename Filter::computations, Expression>;
};

template
<
    typename LaterFilters,
    typename LaterComputations
>
struct reusable_q
{
    template <typename Expression> using fn =
        boost::mp11::mp_or
            <
                boost::mp11::mp_any_of_q
                    <
                        LaterFilters,
                        wanted_by_q<Expression>
                    >,
                boost::mp11::mp_any_of_q
                    <
                        LaterComputations,
                        any_parent_of_q<Expression>
                    >
            >;
};

template
<
    typename AllFilters,
    typename StagedComputations,
    typename Reusables = boost::mp11::mp_list<>,
    typename RemainingFilters = AllFilters,
    typename RemainingStagedComputations =
        boost::mp11::mp_pop_front<StagedComputations>
>
struct all_reusable
{
    using filter = boost::mp11::mp_front<RemainingFilters>;
    using staged_computations =
        boost::mp11::mp_front<RemainingStagedComputations>;
    using later_filters = boost::mp11::mp_pop_front<RemainingFilters>;
    using later_computations =
        boost::mp11::mp_pop_front<RemainingStagedComputations>;
    using reusable = boost::mp11::mp_copy_if_q
        <
            staged_computations,
            reusable_q<later_filters, later_computations>
        >;
    using new_reusables = boost::mp11::mp_push_back<Reusables, reusable>;
public:
    using type = typename all_reusable
        <
            AllFilters,
            StagedComputations,
            new_reusables,
            later_filters,
            later_computations
        >::type;
};

template
<
    typename AllFilters,
    typename StagedComputations,
    typename Reusables
>
struct all_reusable
    <
        AllFilters,
        StagedComputations,
        Reusables,
        boost::mp11::mp_list<>,
        boost::mp11::mp_list<>
    >
{
    using type = Reusables;
};

template
<
    typename Needle,
    typename Haystacks
>
struct index_pair_impl
{
private:
    template <typename List>
    using contains = typename boost::mp11::mp_bind_back
        <
            boost::mp11::mp_contains,
            Needle
        >::fn<List>;
    using outer = boost::mp11::mp_find_if<Haystacks, contains>;
    using inner_list = boost::mp11::mp_at<Haystacks, outer>;
    using inner = boost::mp11::mp_find<inner_list, Needle>;
public:
    using type = std::pair<outer, inner>;
};

template
<
    typename Needle,
    typename Haystacks
>
using index_pair = typename index_pair_impl<Needle, Haystacks>::type;

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_RESULT_PROPAGATION_HPP
