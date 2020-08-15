// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGED_PREDICATE_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGED_PREDICATE_HPP

#include <type_traits>
#include <array>
#include <tuple>

#include <boost/mp11/list.hpp>
#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/set.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/result_propagation.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

struct fail_stage
{
    static constexpr bool stateful = false;

    template
    <
        typename InputList,
        typename StatefulStages,
        typename RemainingStages,
        typename RemainingReusables,
        typename RemainingComputations,
        typename ...InputArr
    >
    static constexpr int staged_apply(const StatefulStages&, const InputArr&...)
    {
        return sign_uncertain;
    }
};

template <typename Stage> using is_stateful = boost::mp11::mp_bool<Stage::stateful>;

template <typename Stage, bool updates = Stage::updates>
struct construct_stage_impl
{
    template <typename Array> static inline Stage apply(const Array& a)
    {
        Stage out(a);
        return out;
    }
};

template <typename Stage>
struct construct_stage_impl<Stage, true>
{
    template <typename Array> static inline Stage apply(const Array&)
    {
        Stage out;
        return out;
    }
};

template <typename Stages>
struct construct_stages_impl
{
    template <typename Array>
    static inline boost::mp11::mp_rename<Stages, std::tuple>
    apply(const Array& a)
    {
        using stage = boost::mp11::mp_front<Stages>;
        std::tuple<stage> first(construct_stage_impl<stage>::apply(a));
        return std::tuple_cat(first,
                              construct_stages_impl
                                    <
                                        boost::mp11::mp_pop_front<Stages>
                                    >::apply(a));
    }
};

template <>
struct construct_stages_impl<std::tuple<>>
{
    template <typename Array>
    static inline std::tuple<> apply(const Array&)
    {
        return std::tuple<>{};
    }
};

template <typename Stage>
using is_updatable = boost::mp11::mp_bool<Stage::updates>;

template <typename UpdatableStages>
struct update_all_impl
{
    template <typename Stages, typename ...Reals>
    static void apply(Stages& stages, const Reals&... args)
    {
        using current_stage = boost::mp11::mp_front<UpdatableStages>;
        std::get<current_stage>(stages).update(args...);
        update_all_impl<boost::mp11::mp_pop_front<UpdatableStages>>
            ::apply(stages, args...);
    }
};

template <>
struct update_all_impl<std::tuple<>>
{
    template <typename Stages, typename ...Reals>
    static void apply(Stages&, const Reals&...) {}
};

template <typename Stage>
using get_arg_count = boost::mp11::mp_size_t<Stage::arg_count>;

template
<
    typename CalculationType,
    typename ...Stages
>
struct staged_predicate
{
public: //TODO: make private
    using ct = CalculationType;
    using stages = std::tuple<Stages...>;
    using stateful_stages = boost::mp11::mp_copy_if<stages, is_stateful>;
    using updatable_stages = boost::mp11::mp_copy_if<stages, is_updatable>;
    stateful_stages m_stages;
    using arg_counts = boost::mp11::mp_transform
        <
            get_arg_count,
            boost::mp11::mp_list<Stages...>
        >;
public:
    static constexpr bool stateful =
        boost::mp11::mp_any_of<stages, is_stateful>::value;
    static constexpr bool updates =
        boost::mp11::mp_any_of<stages, is_updatable>::value;
    static constexpr std::size_t arg_count = boost::mp11::mp_max_element
        <
            arg_counts,
            boost::mp11::mp_less
        >::value;
public: //TODO make private
    using input = argument_list<arg_count>;
    using prev_comps = boost::mp11::mp_list<input>;
    using all_stages = boost::mp11::mp_list<Stages...>;
    using forward_computations =
        typename computations_forward<all_stages, prev_comps>::type;
    using all_re =
        typename all_reusable<all_stages, forward_computations>::type;
public:
    template <typename ...Reals>
    inline staged_predicate(const Reals&... args) : m_stages(
            construct_stages_impl<stateful_stages>::apply(
                std::array<ct, sizeof...(args)>{ static_cast<ct>(args)... }
            )) {}

    template <typename ...Reals>
    inline void update(const Reals&... args)
    {
        update_all_impl<updatable_stages>::apply(m_stages, args...);
    }

    template <typename ...Reals>
    inline int apply(const Reals&... args) const
    {
        std::array<ct, arg_count> input { static_cast<ct>(args)... };
        return next_stage
            <
                prev_comps,
                stateful_stages,
                boost::mp11::mp_push_back<all_stages, fail_stage>,
                all_re,
                boost::mp11::mp_pop_front<forward_computations>
            >::apply(m_stages, input);
    }

/*
    template <typename ...Reals, typename T = boost::mp11::mp_bool<!stateful>>
    static inline typename std::enable_if<T::value, int>::type
        apply(const Reals&... args)
    {
        std::array<ct, arg_count> input { static_cast<ct>(args)... };
        return next_stage
            <
                prev_comps,
                stateful_stages,
                boost::mp11::mp_push_back<all_stages, fail_stage>,
                all_re,
                forward_computations
            >::apply(std::tuple<>{}, input);
    }*/
};

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGED_PREDICATE_HPP
