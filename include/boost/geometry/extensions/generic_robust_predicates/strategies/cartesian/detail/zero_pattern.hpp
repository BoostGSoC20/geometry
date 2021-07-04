#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_ZERO_PATTERN_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_ZERO_PATTERN_HPP

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template <typename Expression>
constexpr int certify_zero(auto input)
{
    if constexpr(Expression::is_leaf)
    {
        if constexpr(Expression::argn != 0)
        {
            return input[Expression::argn - 1] != 0;
        }
        else
        {
            return Expression::value != 0;
        }
    }
    else
    {
        if constexpr(    Expression::operator_type == operator_types::sum
                      || Expression::operator_type == operator_types::difference )
        {
            if constexpr( Expression::left::is_leaf && Expression::right::is_leaf )
            {
                if constexpr( Expression::left::argn != 0 && Expression::right::argn != 0)
                {
                    if constexpr( Expression::operator_type == operator_types::sum )
                    {
                        return   input[ Expression::left::argn - 1 ]
                               + input[ Expression::right::argn - 1 ] != 0;
                    }
                    else
                    {
                        return   input[ Expression::left::argn - 1 ]
                               - input[ Expression::right::argn - 1 ] != 0;
                    }
                }
            }
            else
            {
                return   certify_zero<typename Expression::left>(input)
                       | certify_zero<typename Expression::right>(input);
            }
        }
        if constexpr( Expression::operator_type == operator_types::product )
        {
            return   certify_zero<typename Expression::left>(input)
                   & certify_zero<typename Expression::right>(input);
        }
    }
}

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_ZERO_PATTERN_HPP
