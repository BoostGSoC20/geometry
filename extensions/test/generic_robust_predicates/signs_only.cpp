// Boost.Geometry (aka GGL, Generic Geometry Library)
// Unit Test

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <geometry_test_common.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expressions.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/signs_only_filter.hpp>

using namespace boost::geometry::detail::generic_robust_predicates;

template <typename CalculationType>
void test_all()
{
    using ct = CalculationType;
    using expression = orient2d;
    using filter = signs_only_filter<expression, ct>;
    BOOST_CHECK_EQUAL(-1, filter::apply(2, 1, 0, 1, 1, 2));
    BOOST_CHECK_EQUAL(1, filter::apply(2, 1, 1, 2, 0, 1));
    BOOST_CHECK_EQUAL(sign_uncertain, filter::apply(2, 1, 2, 1, 1, 2));
}

int test_main(int, char* [])
{
    test_all<double>();
    return 0;
}
