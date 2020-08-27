// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SIGNS_ONLY_FILTER_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SIGNS_ONLY_FILTER_HPP

#include <array>

#include <boost/mp11/list.hpp>
#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/set.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/approximate.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template <typename Expression> using is_sign_exact =
    boost::mp11::mp_bool<Expression::sign_exact>;

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename SignArr,
    operator_types Op,
    bool LeftExact,
    bool RightExact,
    typename ...InputArr
>
struct deduce_sign_binary_impl;

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename SignArr,
    typename ...InputArr
>
struct deduce_sign_binary_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        SignArr,
        operator_types::product,
        true,
        true,
        InputArr...
    >
{
    static inline void apply(SignArr& signs,
                             const InputArr&... inputs)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const auto& l = get_approx
                <
                    typename node::left,
                    InputList,
                    Real
                >(inputs...);
        const auto& r = get_approx
                <
                    typename node::right,
                    InputList,
                    Real
                >(inputs...);
        out =   ( l > 0 ? 1 : ( l < 0 ? -1 : 0 ) )
              * ( r > 0 ? 1 : ( r < 0 ? -1 : 0 ) );
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename SignArr,
    typename ...InputArr
>
struct deduce_sign_binary_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        SignArr,
        operator_types::product,
        true,
        false,
        InputArr...
    >
{
    static inline void apply(SignArr& signs,
                             const InputArr&... inputs)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const auto& l = get_approx
                <
                    typename node::left,
                    InputList,
                    Real
                >(inputs...);
        const int& sr =
            signs[boost::mp11::mp_find<All, typename node::right>::value];
        if ( sr == sign_uncertain && l != 0 )
        {
            out = sign_uncertain;
        }
        else
        {
            out = ( l > 0 ? 1 : ( l < 0 ? -1 : 0 ) ) * sr;
        }
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename SignArr,
    typename ...InputArr
>
struct deduce_sign_binary_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        SignArr,
        operator_types::product,
        false,
        true,
        InputArr...
    >
{
    static inline void apply(SignArr& signs,
                             const InputArr&... inputs)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const auto& r = get_approx
                <
                    typename node::right,
                    InputList,
                    Real
                >(inputs...);
        const int& sl =
            signs[boost::mp11::mp_find<All, typename node::left>::value];
        if ( sl == sign_uncertain && r != 0 )
        {
            out = sign_uncertain;
        }
        else
        {
            out = sl * ( r > 0 ? 1 : ( r < 0 ? -1 : 0 ) );
        }
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename SignArr,
    typename ...InputArr
>
struct deduce_sign_binary_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        SignArr,
        operator_types::product,
        false,
        false,
        InputArr...
    >
{
    static inline void apply(SignArr& signs,
                             const InputArr&...)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const int& sl =
            signs[boost::mp11::mp_find<All, typename node::left>::value];
        const int& sr =
            signs[boost::mp11::mp_find<All, typename node::right>::value];
        if ( sl == 0 || sr == 0 )
        {
            out = 0;
        }
        else if ( sl == sign_uncertain || sr == sign_uncertain )
        {
            out = sign_uncertain;
        }
        else
        {
            out = sl * sr;
        }
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename SignArr,
    typename ...InputArr
>
struct deduce_sign_binary_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        SignArr,
        operator_types::sum,
        true,
        true,
        InputArr...
    >
{
    static inline void apply(SignArr& signs,
                             const InputArr&... inputs)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const auto& l = get_approx
                <
                    typename node::left,
                    InputList,
                    Real
                >(inputs...);
        const auto& r = get_approx
                <
                    typename node::right,
                    InputList,
                    Real
                >(inputs...);
        if ( (l > 0 && r >= 0) || (l >= 0 && r > 0) )
        {
            out = 1;
        }
        else if ( (l < 0 && r <= 0) || (l <= 0 && r < 0) )
        {
            out = -1;
        }
        else if ( l == 0 && r == 0 )
        {
            out = 0;
        }
        else
        {
            out = sign_uncertain;
        }
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename SignArr,
    typename ...InputArr
>
struct deduce_sign_binary_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        SignArr,
        operator_types::sum,
        true,
        false,
        InputArr...
    >
{
    static inline void apply(SignArr& signs,
                             const InputArr&... inputs)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const auto& l = get_approx
                <
                    typename node::left,
                    InputList,
                    Real
                >(inputs...);
        const int& sr =
            signs[boost::mp11::mp_find<All, typename node::right>::value];
        if ( sr == sign_uncertain )
        {
            out = sign_uncertain;
        }
        else if ( sr == 1 && l >= 0 )
        {
            out = 1;
        }
        else if ( sr == -1 && l <= 0 )
        {
            out = -1;
        }
        else if ( sr == 0 && l == 0 )
        {
            out = 0;
        }
        else
        {
            out = sign_uncertain;
        }
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename SignArr,
    typename ...InputArr
>
struct deduce_sign_binary_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        SignArr,
        operator_types::sum,
        false,
        true,
        InputArr...
    >
{
    static inline void apply(SignArr& signs,
                             const InputArr&... inputs)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const auto& r = get_approx
                <
                    typename node::right,
                    InputList,
                    Real
                >(inputs...);
        const int& sl =
            signs[boost::mp11::mp_find<All, typename node::left>::value];
        if ( sl == sign_uncertain )
        {
            out = sign_uncertain;
        }
        else if ( sl == 1 && r >= 0 )
        {
            out = 1;
        }
        else if ( sl == -1 && r <= 0 )
        {
            out = -1;
        }
        else if ( sl == 0 && r == 0 )
        {
            out = 0;
        }
        else
        {
            out = sign_uncertain;
        }
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename SignArr,
    typename ...InputArr
>
struct deduce_sign_binary_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        SignArr,
        operator_types::sum,
        false,
        false,
        InputArr...
    >
{
    static inline void apply(SignArr& signs,
                             const InputArr&...)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const int& sl =
            signs[boost::mp11::mp_find<All, typename node::left>::value];
        const int& sr =
            signs[boost::mp11::mp_find<All, typename node::right>::value];
        if ( sl == 0 && sr == 0)
        {
            out = 0;
        }
        else if ( sl == sign_uncertain || sr == sign_uncertain || sl == -sr )
        {
            out = sign_uncertain;
        }
        else if ( sl == sr )
        {
            out = sl;
        }
        else
        {
            out = sl + sr;
        }
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename SignArr,
    typename ...InputArr
>
struct deduce_sign_binary_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        SignArr,
        operator_types::difference,
        true,
        true,
        InputArr...
    >
{
    static inline void apply(SignArr& signs,
                             const InputArr&... inputs)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const auto& l = get_approx
                <
                    typename node::left,
                    InputList,
                    Real
                >(inputs...);
        const auto& r = get_approx
                <
                    typename node::right,
                    InputList,
                    Real
                >(inputs...);
        if ( (l > 0 && r <= 0) || (l >= 0 && r < 0) )
        {
            out = 1;
        }
        else if ( (l < 0 && r >= 0) || (l <= 0 && r > 0) )
        {
            out = -1;
        }
        else if ( l == 0 && r == 0 )
        {
            out = 0;
        }
        else
        {
            out = sign_uncertain;
        }
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename SignArr,
    typename ...InputArr
>
struct deduce_sign_binary_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        SignArr,
        operator_types::difference,
        true,
        false,
        InputArr...
    >
{
    static inline void apply(SignArr& signs,
                             const InputArr&... inputs)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const auto& l = get_approx
                <
                    typename node::left,
                    InputList,
                    Real
                >(inputs...);
        const int& sr =
            signs[boost::mp11::mp_find<All, typename node::right>::value];
        if ( sr == sign_uncertain )
        {
            out = sign_uncertain;
        }
        else if ( sr == -1 && l >= 0 )
        {
            out = 1;
        }
        else if ( sr == 1 && l <= 0 )
        {
            out = -1;
        }
        else if ( sr == 0 && l == 0 )
        {
            out = 0;
        }
        else
        {
            out = sign_uncertain;
        }
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename SignArr,
    typename ...InputArr
>
struct deduce_sign_binary_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        SignArr,
        operator_types::difference,
        false,
        true,
        InputArr...
    >
{
    static inline void apply(SignArr& signs,
                             const InputArr&... inputs)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const auto& r = get_approx
                <
                    typename node::right,
                    InputList,
                    Real
                >(inputs...);
        const int& sl =
            signs[boost::mp11::mp_find<All, typename node::left>::value];
        if ( sl == sign_uncertain )
        {
            out = sign_uncertain;
        }
        else if ( sl == 1 && r <= 0 )
        {
            out = 1;
        }
        else if ( sl == -1 && r >= 0 )
        {
            out = -1;
        }
        else if ( sl == 0 && r == 0 )
        {
            out = 0;
        }
        else
        {
            out = sign_uncertain;
        }
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename SignArr,
    typename ...InputArr
>
struct deduce_sign_binary_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        SignArr,
        operator_types::difference,
        false,
        false,
        InputArr...
    >
{
    static inline void apply(SignArr& signs,
                             const InputArr&...)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const int& sl =
            signs[boost::mp11::mp_find<All, typename node::left>::value];
        const int& sr =
            signs[boost::mp11::mp_find<All, typename node::right>::value];
        if ( sl == 0 && sr == 0 )
        {
            out = 0;
        }
        else if ( sl == sr || sl == sign_uncertain || sr == sign_uncertain )
        {
            out = sign_uncertain;
        }
        else if ( sl == 0 || sr == 0 )
        {
            out = sl - sr;
        }
        else
        {
            out = sl;
        }
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename SignArr,
    operator_arities Arity,
    typename ...InputArr
>
struct deduce_sign_arity_helper_impl {};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename SignArr,
    typename ...InputArr
>
struct deduce_sign_arity_helper_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        SignArr,
        operator_arities::binary,
        InputArr...
    >
{
    static inline void apply(SignArr& s, const InputArr&... inputs)
    {
        using node = boost::mp11::mp_front<Remaining>;
        deduce_sign_binary_impl
            <
                All,
                Remaining,
                InputList,
                Real,
                SignArr,
                node::operator_type,
                node::left::sign_exact,
                node::right::sign_exact,
                InputArr...
            >::apply(s, inputs...);
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename SignArr,
    std::size_t RemainingSize,
    typename ...InputArr
>
struct deduce_signs_remainder_impl
{
    static inline void apply(SignArr& signs,
                             const InputArr&... inputs)
    {
        using node = boost::mp11::mp_front<Remaining>;
        deduce_sign_arity_helper_impl
            <
                All,
                Remaining,
                InputList,
                Real,
                SignArr,
                node::operator_arity,
                InputArr...
            >::apply(signs, inputs...);
        deduce_signs_remainder_impl
            <
                All,
                boost::mp11::mp_pop_front<Remaining>,
                InputList,
                Real,
                SignArr,
                boost::mp11::mp_size<Remaining>::value - 1,
                InputArr...
            >::apply(signs, inputs...);
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename SignArr,
    typename ...InputArr
>
struct deduce_signs_remainder_impl
    <
        All,
        Remaining,
        InputList,
        Real,
        SignArr,
        1,
        InputArr...
    >
{
    static inline void apply(SignArr& signs,
                             const InputArr&... inputs)
    {
        using node = boost::mp11::mp_front<Remaining>;
        deduce_sign_arity_helper_impl
            <
                All,
                Remaining,
                InputList,
                Real,
                SignArr,
                node::operator_arity,
                InputArr...
            >::apply(signs, inputs...);
    }
};

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename SignArr,
    typename ...InputArr
>
inline void deduce_signs_remainder(SignArr& signs,
                                  const InputArr&... inputs)
{
    deduce_signs_remainder_impl
        <
            All,
            Remaining,
            InputList,
            Real,
            SignArr,
            boost::mp11::mp_size<Remaining>::value,
            InputArr...
        >::apply(signs, inputs...);
}

template
<
    typename All,
    typename Remaining,
    typename InputList,
    typename Real,
    typename SignArr,
    typename ...InputArr
>
inline void deduce_signs(SignArr& signs,
                         const InputArr&... inputs)
{
    deduce_signs_remainder
        <
            All,
            Remaining,
            InputList,
            Real
        >(signs, inputs...);
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
        using arg_list_input = argument_list<sizeof...(Reals)>;
        using arg_list = boost::mp11::mp_list<evals_sign_exact, arg_list_input>;
        std::array<Real, sizeof...(args)> input {{ static_cast<Real>(args)... }};
        std::array<Real, boost::mp11::mp_size<evals_sign_exact>::value>
            results_sign_exact;
        approximate_interim
            <
                evals_sign_exact,
                arg_list,
                Real
            >(results_sign_exact, input);
        using allm = boost::mp11::mp_push_front<arg_list, evals_sign_exact>;
        std::array<int, boost::mp11::mp_size<non_exact_signs>::value>
            remainder_signs;
        deduce_signs
            <
                non_exact_signs,
                non_exact_signs,
                allm,
                Real
            >(remainder_signs, results_sign_exact, input);
        return remainder_signs[
            boost::mp11::mp_find<non_exact_signs, Expression>::value];
    }
};

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SIGNS_ONLY_FILTER_HPP
