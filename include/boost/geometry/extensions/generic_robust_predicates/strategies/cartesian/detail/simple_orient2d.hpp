// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SIMPLE_ORIENT2D_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SIMPLE_ORIENT2D_HPP

#include <cstdint>
#include <limits>
#include <cmath>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expressions.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/almost_static_filter.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

//Once std math functions become constexpr, this could be made entirely generic
template <int digits>
struct phi_impl
{};

template <> struct phi_impl<53>
{
    static constexpr int_least64_t value = 94'906'264;
};

template <> struct phi_impl<24>
{
    static constexpr int_least32_t value = 4'094;
};

template <> struct phi_impl<64>
{
    static constexpr int_least64_t value = 4'294'967'294;
};

template <> struct phi_impl<11>
{
    static constexpr int_least16_t value = 44;
};

template <typename CalculationType>
static constexpr auto phi =
    phi_impl<std::numeric_limits<CalculationType>::digits>::value;

template <typename CalculationType>
static constexpr CalculationType theta =
      3 * std::numeric_limits<CalculationType>::epsilon() / 2
    - (phi<CalculationType> - 22) * (std::numeric_limits<CalculationType>::epsilon() / 2)
                                  * (std::numeric_limits<CalculationType>::epsilon() / 2);

template <typename CalculationType>
inline CalculationType ufp(CalculationType a)
{
    if(a == 0.)
    {
        return 0.;
    }
    else
    {
        return std::exp2(std::floor(std::log2(std::abs(a))));
    }
}

// The following filters are based on "Simple Floating-Point Filters for the
// Two-Dimensional Orientation Problem" by Ozaki et al.

template <typename CalculationType = double>
struct simple_orient2d_semi_static
{
private:
    using ct = CalculationType;
public:
    static constexpr bool stateful = false;
    static constexpr bool updates = false;

    template <typename ...Reals>
    inline simple_orient2d_semi_static(const Reals&...) {}

    template <typename CT>
    static inline ct error_bound(CT ax, CT ay, CT bx, CT by, CT cx, CT cy)
    {
        const ct l = (ax - cx) * (by - cy);
        const ct r = (bx - cx) * (ay - cy);
        return   theta<ct>
               * ( std::abs(l + r) + std::numeric_limits<ct>::min() );
    }

    template <typename CT>
    static inline int apply(CT ax, CT ay, CT bx, CT by, CT cx, CT cy)
    {
        const ct l = (ax - cx) * (by - cy);
        const ct r = (bx - cx) * (ay - cy);
        const ct det = l - r;
        const ct errbound =
              theta<ct>
            * ( std::abs(l + r) + std::numeric_limits<ct>::min() );
        if( det > errbound)
        {
            return 1;
        }
        else if( det < -errbound )
        {
            return -1;
        }
        else if ( det == 0 && errbound == 0 )
        {
            return 0;
        }
        return sign_uncertain;
    }
};

template <typename CalculationType>
struct simple_orient2d_static
{
private:
    using ct = CalculationType;
    ct m_error_bound;
public:
    static constexpr bool stateful = true;
    static constexpr bool updates = false;

    inline simple_orient2d_static() : m_error_bound(0) {}

    inline ct error_bound() const { return m_error_bound; }

    static inline ct compute_eb(const std::array<ct, 12>& a)
    {
        const ct mx = std::max(std::max(a[0], a[2]), a[4]);
        const ct my = std::max(std::max(a[1], a[3]), a[5]);
        const ct nx = std::min(std::min(a[6], a[8]), a[10]);
        const ct ny = std::min(std::min(a[7], a[9]), a[11]);
        const ct alpha = mx - nx;
        const ct beta  = my - ny;
        constexpr ct u = std::numeric_limits<ct>::epsilon() / 2;
        const ct T_2 =   2 * alpha * u * ufp(beta) + 2 * beta * u * ufp(alpha)
                       + 2 * u * ufp( alpha * beta )
                       + 2 * u * u * ufp(alpha) * ufp(beta);
        return std::nextafter(T_2 + 3 * u * ufp(T_2),
                              std::numeric_limits<ct>::max());
    }

    template <typename ...Reals>
    inline simple_orient2d_static(const Reals&... args)
        : m_error_bound(compute_eb(std::array<ct, 12>{static_cast<ct>(args)...}))
    {
        static_assert(sizeof...(Reals) == 12,
                      "Number of constructor arguments is incompatible with error expression.");
    }

    inline simple_orient2d_static(const std::array<ct, 12>& extrema)
        : m_error_bound(compute_eb(extrema)) {}

    inline int apply(ct ax, ct ay, ct bx, ct by, ct cx, ct cy) const
    {
        const ct l = (ax - cx) * (by - cy);
        const ct r = (bx - cx) * (ay - cy);
        const ct det = l - r;
        if(det > m_error_bound)
        {
            return 1;
        }
        else if(det < -m_error_bound)
        {
            return -1;
        }
        else if(m_error_bound == 0 && det == 0)
        {
            return 0;
        }
        return sign_uncertain;
    }
};

template <typename CalculationType>
using simple_orient2d_almost_static = almost_static_filter
    <
        orient<2>,
        CalculationType,
        simple_orient2d_static<CalculationType>
    >;

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SIMPLE_ORIENT2D_HPP
