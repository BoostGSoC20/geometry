// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SIGNS_ONLY_FILTER_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SIGNS_ONLY_FILTER_HPP

#include <type_traits>
#include <array>

#include <boost/mp11/list.hpp>
#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/set.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_eval.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template <typename Expression> using is_sign_exact =
    boost::mp11::mp_bool<Expression::sign_exact>;

template <typename Real>
constexpr std::enable_if_t<std::is_floating_point<Real>::value, int> sign(const Real a)
{
    if( a > 0 )
    {
        return 1;
    }
    else if ( a < 0 )
    {
        return -1;
    }
    else
    {
        return 0;
    }
}

template <operator_types Op>
constexpr int sign(const int a, const int b);

template <>
constexpr int sign<operator_types::sum>(const int a, const int b)
{
    if(a == 0)
    {
        return b;
    }
    else if(b == 0)
    {
        return a;
    }
    else if(a == b && a != sign_uncertain)
    {
        return a;
    }
    else
    {
        return sign_uncertain;
    }
}

template <>
constexpr int sign<operator_types::difference>(const int a, const int b)
{
    if(b == 0)
    {
        return a;
    }
    else if(b == -a)
    {
        return a;
    }
    else if(a == 0 && b != sign_uncertain)
    {
        return -b;
    }
    else
    {
        return sign_uncertain;
    }
}

template <>
constexpr int sign<operator_types::product>(const int a, const int b)
{
    if(a == 0 || b == 0)
    {
        return 0;
    }
    else if(a != sign_uncertain && b != sign_uncertain)
    {
        return a * b;
    }
    else
    {
        return sign_uncertain;
    }
}

template
<
    typename Expression,
    typename ExactSignExpressions,
    typename DeducedSignExpressions,
    bool IsExact = Expression::sign_exact,
    bool IsLeaf = Expression::is_leaf
>
struct get_sign {};

template
<
    typename Expression,
    typename ExactSignExpressions,
    typename DeducedSignExpressions
>
struct get_sign
    <
        Expression,
        ExactSignExpressions,
        DeducedSignExpressions,
        false,
        false
    >
{
    template <typename InputArr, typename ExactArr, typename SignArr>
    static constexpr int apply(const InputArr&, const ExactArr&, const SignArr& s)
    {
        return s[boost::mp11::mp_find<DeducedSignExpressions, Expression>::value];
    }
};

template
<
    typename Expression,
    typename ExactSignExpressions,
    typename DeducedSignExpressions
>
struct get_sign
    <
        Expression,
        ExactSignExpressions,
        DeducedSignExpressions,
        true,
        false
    >   
{
    template <typename InputArr, typename ExactArr, typename SignArr>
    static constexpr int apply(const InputArr&, const ExactArr& e, const SignArr&)
    {
        return sign(e[boost::mp11::mp_find<ExactSignExpressions, Expression>::value]);
    }
};

template
<
    typename Expression,
    typename ExactSignExpressions,
    typename DeducedSignExpressions
>
struct get_sign
    <
        Expression,
        ExactSignExpressions,
        DeducedSignExpressions,
        true,
        true
    >
{
    template <typename InputArr, typename ExactArr, typename SignArr>
    static constexpr int apply(const InputArr& input, const ExactArr&, const SignArr&)
    {
        return sign(get_arg_or_const<Expression>(input));
    }
};

template
<
    typename Expression,
    typename ExactSignExpressions,
    typename DeducedSignExpressions
>
struct deduce_sign_binary
{
private:
    using left = typename Expression::left;
    using right = typename Expression::right;
public:
    template <typename InputArr, typename ExactArr, typename SignArr>
    static constexpr void apply(const InputArr& i,
                                const ExactArr& e,
                                SignArr& s)
    {
        s[boost::mp11::mp_find<DeducedSignExpressions, Expression>::value] =
            sign<Expression::operator_type>(
                    get_sign
                        <
                            left,
                            ExactSignExpressions,
                            DeducedSignExpressions
                        >::apply(i, e, s),
                    get_sign
                        <
                            right,
                            ExactSignExpressions,
                            DeducedSignExpressions
                        >::apply(i, e, s));
    }
};

template
<
    typename ExactSignExpressions,
    typename InputArr,
    typename ExactArr,
    typename SignArr,
    typename ...DeducedSignExpressions
>
constexpr void deduce_signs(const InputArr& i,
                            const ExactArr& e,
                            SignArr& s,
                            boost::mp11::mp_list<DeducedSignExpressions...>)
{
    using deduced_list = boost::mp11::mp_list<DeducedSignExpressions...>;
    using dummy = int[];
    (void)dummy{
        0,
        (deduce_sign_binary
            <
                DeducedSignExpressions,
                ExactSignExpressions,
                deduced_list
            >::apply(i, e, s), 0)...
    };
}

// The following filter tries to deduce the sign of an expression solely based
// on the signs of its subexpressions, e.g. we know that 
// subexpression1 - subexpression2 > 0 if subexpression1 > 0 and
// subexpression2 <= 0 or subexpression1 >= 0 and subexpression2 < 0. It makes
// use of floating point approximation for all subexpressions for which the
// approximations are guaranteed to have the correct sign, e.g. a - b and a + b
// have the correct sign in floating point arithmetic, if a and b are exact
// input values and subexpression1 * subexpression2 is guaranteed to have the
// correct sign in floating point arithmetic, if the same holds for each
// subexpression.

template <typename Expression, typename Real>
struct signs_only_filter
{
private:
    using non_exact_signs_po =
        typename boost::mp11::mp_unique<post_order<Expression, is_sign_exact>>;
    using non_exact_signs =
        typename boost::mp11::mp_remove_if<non_exact_signs_po, is_sign_exact>;
    using stack = typename boost::mp11::mp_unique<post_order<Expression>>;
    using evals = typename boost::mp11::mp_remove_if<stack, is_leaf>;
    using evals_sign_exact =
        typename boost::mp11::mp_copy_if<evals, is_sign_exact>;
    using ct = Real;
public:
    static constexpr bool stateful = false;
    static constexpr bool updates = false;
    static constexpr std::size_t arg_count = max_argn<Expression>::value;

    template <typename ...Reals>
    static inline int apply(const Reals&... args)
    {
        std::array<ct, sizeof...(args)> input {{ static_cast<Real>(args)... }};
        std::array<ct, boost::mp11::mp_size<evals_sign_exact>::value>
            results_sign_exact;
        evaluate_expressions(input, results_sign_exact, evals_sign_exact{});
        std::array<int, boost::mp11::mp_size<non_exact_signs>::value>
            remainder_signs;
        deduce_signs
            <
                evals_sign_exact
            >(input, results_sign_exact, remainder_signs, non_exact_signs{});
        return get_sign
            <
                Expression,
                evals_sign_exact,
                non_exact_signs
            >::apply(input, results_sign_exact, remainder_signs);
    }
};

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SIGNS_ONLY_FILTER_HPP
