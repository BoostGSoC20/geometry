// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_PREDICATE_APPROXIMATION_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_PREDICATE_APPROXIMATION_HPP

#include <array>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_eval.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template
<
    typename Expression,
    typename CT,
    operator_types Op = Expression::operator_type
>
struct approx_sign
{
    template <typename ...Reals>
    static inline int apply(const Reals&... args)
    {
        std::array<CT, sizeof...(Reals)> input {static_cast<CT>(args)...};
        const CT approx = evaluate_expression<Expression>(input);
        if(approx > 0)
        {
            return 1;
        }
        else if(approx < 0)
        {
            return -1;
        }
        else
        {
            return 0;
        }
    }
};

template <typename Expression, typename CT>
struct approx_sign<Expression, CT, operator_types::difference>
{
    template <typename ...Reals>
    static inline int apply(const Reals&... args)
    {
        std::array<CT, sizeof...(Reals)> input {static_cast<CT>(args)...};
        const CT left_approx =
            evaluate_expression<typename Expression::left>(input);
        const CT right_approx =
            evaluate_expression<typename Expression::right>(input);
        if(left_approx > right_approx)
        {
            return 1;
        }
        else if(left_approx < right_approx)
        {
            return -1;
        }
        else
        {
            return 0;
        }
    }
};

template <typename Expression, typename CT>
struct approx_sign<Expression, CT, operator_types::sum>
{
    template <typename ...Reals>
    static inline int apply(const Reals&... args)
    {
        std::array<CT, sizeof...(Reals)> input {static_cast<CT>(args)...};
        const CT left_approx =
            evaluate_expression<typename Expression::left>(input);
        const CT right_approx =
            -evaluate_expression<typename Expression::right>(input);
        if(left_approx > right_approx)
        {
            return 1;
        }
        else if(left_approx < right_approx)
        {
            return -1;
        }
        else
        {
            return 0;
        }
    }
};

template
<
    typename Expression,
    typename CalculationType
>
struct predicate_approximation
{
    static constexpr bool stateful = false;
    static constexpr bool updates = false;

    template <typename ...Reals>
    static inline int apply(const Reals&... args)
    {
        return approx_sign<Expression, CalculationType>::apply(args...);
    }

    template <typename ...Reals>
    inline int operator()(Reals const&... args) const
    {
        return apply(args...);
    }
};

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_PREDICATE_APPROXIMATION_HPP
