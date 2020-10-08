// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPANSION_EVAL_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPANSION_EVAL_HPP

#include <cstddef>
#include <array>

#include <boost/mp11/integral.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/algorithm.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expansion_arithmetic.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

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

template <int> using no_zero_elimination_policy = boost::mp11::mp_false;

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
    static constexpr Iter apply(Args...args)
    {
        return expansion_plus
            <
                LeftLength,
                RightLength,
                Inplace,
                no_zero_elimination_policy
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
        operator_types::difference,
        LeftLength,
        RightLength,
        Inplace,
        Iter,
        StageB
    >
{
    template<typename ...Args>
    static constexpr Iter apply(Args...args)
    {
        return expansion_minus
            <
                LeftLength,
                RightLength,
                Inplace,
                StageB,
                no_zero_elimination_policy
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
    static constexpr Iter apply(Args...args)
    {
        return expansion_times
            <
                LeftLength,
                RightLength,
                Inplace,
                no_zero_elimination_policy
            >(args...);
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
    static constexpr Iter apply(Iter begin, Iter, const Reals&... args)
    {
        std::array<Real, sizeof...(Reals)> input
            {{ static_cast<Real>(args)... }};
        Real left_val = input[left::argn - 1];
        Real right_val = input[right::argn - 1];
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
    static constexpr Iter apply(Iter begin, Iter, const Reals&... args)
    {
        std::array<Real, sizeof...(Reals)> input
            {{ static_cast<Real>(args)... }};
        Real left_val = input[left::argn - 1];
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
    static constexpr Iter apply(Iter begin, Iter, const Reals&... args)
    {
        std::array<Real, sizeof...(Reals)> input
            {{ static_cast<Real>(args)... }};
        Real right_val = input[right::argn - 1];
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
    static constexpr Iter apply(Iter begin, Iter, const Reals&...)
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
    static constexpr Iter apply(Iter begin, Iter end, const Reals&... args)
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
    static constexpr Iter apply(Iter begin, Iter end, const Reals&... args)
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

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPANSION_EVAL_HPP
