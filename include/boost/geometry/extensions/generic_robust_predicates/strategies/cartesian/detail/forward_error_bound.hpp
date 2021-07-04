#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_FORWARD_ERROR_BOUND_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_FORWARD_ERROR_BOUND_HPP

#include <array>
#include <limits>

#include <boost/mp11/utility.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/simple_orient2d.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_a_error_bound.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template <typename Expression, typename CalculationType, typename Rules>
struct forward_error_bound;


template <typename CalculationType>
struct underflow_guard_constant : public static_constant_interface<CalculationType>
{
    // 2 * u_N
    static constexpr CalculationType value =
        2 * std::numeric_limits<CalculationType>::min();
    static constexpr bool non_negative = true;
};

template
<
    bool UnderflowProtection = true
>
struct ozaki_simple_fp_lemma_31
{
    template <typename Expression, typename>
    static constexpr bool applicable()
    {
        if constexpr( Expression::operator_type == operator_types::product )
        {
            if constexpr(    (   Expression::left::operator_type == operator_types::sum
                              || Expression::left::operator_type == operator_types::difference)
                          && (   Expression::right::operator_type == operator_types::sum
                              || Expression::right::operator_type == operator_types::difference) )
            {
                return Expression::left::left::is_leaf && Expression::left::right::is_leaf
                    && Expression::right::left::is_leaf && Expression::right::right::is_leaf;
            }
        }
        return false;
    }

    template <typename Expression, typename CalculationType, typename>
    struct error_bound
    {
        using magnitude = 
            mp11::mp_if_c
                <
                    UnderflowProtection,
                    sum
                        <
                            abs<Expression>,
                            underflow_guard_constant<CalculationType>
                        >,
                    abs<Expression>
                >;

        static constexpr std::array<int, 3> a {3, -(phi<CalculationType> - 14), 0};
    };
};

struct exact_leaves
{
    template <typename Expression, typename>
    static constexpr bool applicable()
    {
        return Expression::is_leaf;
    }

    template <typename Expression, typename, typename>
    struct error_bound
    {
        using magnitude = abs<Expression>;
        static constexpr std::array<int, 3> a {0, 0, 0};
    };
};

struct exact_leaves_sumdiff
{
    template <typename Expression, typename>
    static constexpr bool applicable()
    {
        if constexpr(    Expression::operator_type == operator_types::sum 
                      || Expression::operator_type == operator_types::difference )
        {
            return Expression::left::is_leaf && Expression::right::is_leaf;
        }
        return false;
    }

    template <typename Expression, typename, typename>
    struct error_bound
    {
        using magnitude = abs<Expression>;
        static constexpr std::array<int, 3> a {1, 0, 0};
    };
};

template
<
    bool UnderflowProtection = true
>
struct exact_leaves_product
{
    template <typename Expression, typename>
    static constexpr bool applicable()
    {
        if constexpr( Expression::operator_type == operator_types::product )
        {
            return Expression::left::is_leaf && Expression::right::is_leaf;
        }
        return false;
    }

    template <typename Expression, typename CalculationType, typename>
    struct error_bound
    {
        using magnitude =
            mp11::mp_if_c
                <
                    UnderflowProtection,
                    sum
                        <
                            abs<Expression>,
                            underflow_guard_constant<CalculationType>
                        >,
                    abs<Expression>
                >;
        static constexpr std::array<int, 3> a {1, 0, 0};
    };
};

struct inexacts_sumdiff
{
    template <typename Expression, typename>
    static constexpr bool applicable()
    {
        return    Expression::operator_type == operator_types::sum
               || Expression::operator_type == operator_types::difference;
    }

    template <typename Expression, typename CalculationType, typename Rules>
    struct error_bound
    {
    private:
        using leb = typename forward_error_bound<typename Expression::left, CalculationType, Rules>::error_bound;
        using reb = typename forward_error_bound<typename Expression::right, CalculationType, Rules>::error_bound;
        static constexpr std::array<int, 3> max_a = coeff_max(leb::a, reb::a);
    public:
        using magnitude = sum<typename leb::magnitude, typename reb::magnitude>;
        static constexpr std::array<int, 3> a = coeff_inc_first(coeff_mult_by_1_plus_eps(max_a));
    };
};

template
<
    bool UnderflowProtection = true
>
struct inexacts_product
{
    template <typename Expression, typename>
    static constexpr bool applicable()
    {
        return Expression::operator_type == operator_types::product;
    }

    template <typename Expression, typename CalculationType, typename Rules>
    struct error_bound
    {
    private:
        using leb = typename forward_error_bound
            <
                typename Expression::left,
                CalculationType,
                Rules
            >::error_bound;
        using reb = typename forward_error_bound
            <
                typename Expression::right,
                CalculationType,
                Rules
            >::error_bound;
        static constexpr std::array<int, 3> a_prod = coeff_product(leb::a, reb::a);
    public:
        using magnitude = mp11::mp_if_c
            <
                UnderflowProtection,
                sum
                    <
                        product
                            <
                                typename leb::magnitude,
                                typename reb::magnitude
                            >,
                        underflow_guard_constant<CalculationType>
                    >,
                product<typename leb::magnitude, typename reb::magnitude>
            >;
        static constexpr std::array<int, 3> a = coeff_inc_first(coeff_mult_by_1_plus_eps(a_prod));
    };
};

template <typename Expression, typename CalculationType>
struct applicable_to
{
    template <typename Rule>
    using fn = mp11::mp_bool<Rule::template applicable<Expression, CalculationType>()>;
};

template <typename Expression, typename CalculationType, typename Rules>
struct forward_error_bound
{
    using rule = mp11::mp_at
        <
            Rules,
            mp11::mp_find_if_q
                <
                    Rules,
                    applicable_to<Expression, CalculationType>
                >
        >;

    using error_bound = typename rule::template error_bound
        <
            Expression,
            CalculationType,
            Rules
        >;
};

template <typename Expression, typename CalculationType, typename Rules>
struct forward_error_condition_sumdiff
{
private:
    using leb = typename forward_error_bound<typename Expression::left, CalculationType, Rules>::error_bound;
    using reb = typename forward_error_bound<typename Expression::right, CalculationType, Rules>::error_bound;
    static constexpr auto max_a = coeff_max(leb::a, reb::a);
    static constexpr std::array<int, 3> a =
        coeff_mult_by_1_plus_eps(coeff_mult_by_1_plus_eps(coeff_div_by_1_minus_eps(max_a)));
    static constexpr int c = round_to_next_2_pow<a[0]>();
    static constexpr int eps_square_coeff =
        a[2] > 0 ? c * ((a[1] + 1) / c + 1) : c * (a[1] / c + 1);
public:
    using magnitude = sum<typename leb::magnitude, typename reb::magnitude>;
    static constexpr std::array<int, 2> coefficients {a[0], eps_square_coeff};
};

template
<
    typename Expression,
    typename CalculationType,
    typename Rules
>
struct forward_error_bound_expression_impl
{
private:
    using ct = CalculationType;
    using fe_cond = forward_error_condition_sumdiff<Expression, CalculationType, Rules>;
    static constexpr ct eps = std::numeric_limits<ct>::epsilon() / 2.0;
    struct constant : public static_constant_interface<ct>
    {
        static constexpr ct value =
              fe_cond::coefficients[0] * eps
            + fe_cond::coefficients[1] * eps * eps;
        static constexpr bool non_negative = true;
    };
public:
    using type = product<constant, typename fe_cond::magnitude>;
};

template
<
    typename Expression,
    typename CalculationType,
    typename Rules
>
using forward_error_bound_expression =
    typename forward_error_bound_expression_impl<Expression, CalculationType, Rules>::type;

using all_rules = mp11::mp_list<exact_leaves, exact_leaves_sumdiff, exact_leaves_product<false>, inexacts_sumdiff, ozaki_simple_fp_lemma_31<false>, inexacts_product<false>>;
using all_rules_up = mp11::mp_list<exact_leaves, exact_leaves_sumdiff, exact_leaves_product<true>, inexacts_sumdiff, ozaki_simple_fp_lemma_31<true>, inexacts_product<true>>;

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_FORWARD_ERROR_BOUND_HPP
