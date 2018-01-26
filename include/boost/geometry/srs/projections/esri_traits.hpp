// Boost.Geometry

// Copyright (c) 2017, Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_PROJECTIONS_ESRI_TRAITS_HPP
#define BOOST_GEOMETRY_PROJECTIONS_ESRI_TRAITS_HPP


#include <boost/geometry/core/tags.hpp>

#include <boost/geometry/srs/projections/factory.hpp>
#include <boost/geometry/srs/projections/impl/projects.hpp>
#include <boost/geometry/srs/projections/srid_traits.hpp>


namespace boost { namespace geometry { namespace projections
{

#ifndef DOXYGEN_NO_DETAIL
namespace detail
{

/*!
    \brief ESRI traits
    \details With help of the ESRI traits library users can statically use projections
        or coordinate systems specifying an ESRI code. The correct projections for transformations
        are used automically then, still keeping static polymorphism.
    \ingroup projection
    \tparam ESRI esri code
*/
template <size_t ESRI>
struct esri_traits
{
    // Specializations define:
    // - type to get projection type
    // - function par to get parameters
};

BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37001, longlat, WGS66, "+proj=longlat +ellps=WGS66 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 37002, longlat, 6378166, 6356784.283607107, "+proj=longlat +a=6378166 +b=6356784.283607107 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 37003, longlat, 6378150, 6356768.337244385, "+proj=longlat +a=6378150 +b=6356768.337244385 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37004, longlat, fschr60m, "+proj=longlat +ellps=fschr60m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 37005, longlat, 6378270, 6356794.343434343, "+proj=longlat +a=6378270 +b=6356794.343434343 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 37006, longlat, 6377295.664, 6356094.667915204, "+proj=longlat +a=6377295.664 +b=6356094.667915204 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 37007, longlat, 6376896, 6355834.846687363, "+proj=longlat +a=6376896 +b=6355834.846687363 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 37008, longlat, 6370997, "+proj=longlat +a=6370997 +b=6370997 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37201, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 37202, longlat, 6377276.345, 6356075.41314024, "+proj=longlat +a=6377276.345 +b=6356075.41314024 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 37203, longlat, 6377301.243, 6356100.230165384, "+proj=longlat +a=6377301.243 +b=6356100.230165384 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37204, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37205, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37206, longlat, clrk80, "+proj=longlat +ellps=clrk80 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37207, longlat, fschr60m, "+proj=longlat +ellps=fschr60m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37208, longlat, clrk80, "+proj=longlat +ellps=clrk80 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37211, longlat, clrk80, "+proj=longlat +ellps=clrk80 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37212, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37213, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37214, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37215, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37216, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37217, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37218, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37219, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37220, longlat, clrk66, "+proj=longlat +ellps=clrk66 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37221, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37222, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 37223, longlat, 6378249.2, 6356514.999904194, "+proj=longlat +a=6378249.2 +b=6356514.999904194 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37224, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37226, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37227, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37228, longlat, clrk80, "+proj=longlat +ellps=clrk80 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 37229, longlat, 6378270, 6356794.343434343, "+proj=longlat +a=6378270 +b=6356794.343434343 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37230, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37231, longlat, aust_SA, "+proj=longlat +ellps=aust_SA +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37232, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37233, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37234, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37235, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37237, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37238, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37239, longlat, clrk66, "+proj=longlat +ellps=clrk66 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37240, longlat, clrk80, "+proj=longlat +ellps=clrk80 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37241, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37242, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37243, longlat, clrk66, "+proj=longlat +ellps=clrk66 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37245, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37246, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37247, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37249, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37250, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37251, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37252, longlat, clrk66, "+proj=longlat +ellps=clrk66 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37253, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37254, longlat, clrk80, "+proj=longlat +ellps=clrk80 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37255, longlat, bessel, "+proj=longlat +ellps=bessel +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37257, longlat, krass, "+proj=longlat +ellps=krass +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37259, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 37260, longlat, clrk66, "+proj=longlat +ellps=clrk66 +no_defs");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 53001, "");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 53002, eqc, 6371000, "+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 53003, mill, 6371000, "+proj=mill +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +R_A +a=6371000 +b=6371000 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 53004, merc, 6371000, "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 53008, sinu, 6371000, "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 53009, moll, 6371000, "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 53010, eck6, 6371000, "+proj=eck6 +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 53011, "");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 53012, eck4, 6371000, "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 53013, "");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 53014, "");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 53015, "");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 53016, gall, 6371000, "+proj=gall +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 53017, "");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 53018, "");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 53019, "");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 53021, poly, 6371000, "+proj=poly +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 53022, "");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 53023, "");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 53024, bonne, 6371000, "+proj=bonne +lon_0=0 +lat_1=60 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 53025, omerc, 6371000, "+proj=omerc +lat_0=40 +lon_1=0 +lat_1=0 +lon_2=0 +lat_2=0 +k=1 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 53026, stere, 6371000, "+proj=stere +lat_0=0 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 53027, eqdc, 6371000, "+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=60 +lat_2=60 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 53028, cass, 6371000, "+proj=cass +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 53029, vandg, 6371000, "+proj=vandg +lon_0=0 +x_0=0 +y_0=0 +R_A +a=6371000 +b=6371000 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 53030, robin, 6371000, "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 53031, tpeqd, 6371000, "+proj=tpeqd +lat_1=0 +lon_1=0 +lat_2=60 +lon_2=60 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 53032, aeqd, 6371000, "+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 54001, "");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 54002, eqc, WGS84, WGS84, "+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 54003, mill, 6371007.1810824294, "+proj=mill +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +R_A +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 54004, merc, WGS84, WGS84, "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 54008, sinu, WGS84, WGS84, "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 54009, moll, WGS84, WGS84, "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 54010, eck6, WGS84, WGS84, "+proj=eck6 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 54011, "");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 54012, eck4, WGS84, WGS84, "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 54013, "");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 54014, "");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 54015, "");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 54016, gall, WGS84, WGS84, "+proj=gall +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 54017, "");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 54018, "");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 54019, "");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 54021, poly, WGS84, WGS84, "+proj=poly +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 54022, "");
//BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS (esri, 54023, "");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 54024, bonne, WGS84, WGS84, "+proj=bonne +lon_0=0 +lat_1=60 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 54025, omerc, WGS84, WGS84, "+proj=omerc +lat_0=40 +lon_1=0 +lat_1=0 +lon_2=0 +lat_2=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 54026, stere, WGS84, WGS84, "+proj=stere +lat_0=0 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 54027, eqdc, WGS84, WGS84, "+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=60 +lat_2=60 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 54028, cass, WGS84, WGS84, "+proj=cass +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_S (esri, 54029, vandg, 6371007.1810824294, "+proj=vandg +lon_0=0 +x_0=0 +y_0=0 +R_A +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 54030, robin, WGS84, WGS84, "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 54031, tpeqd, WGS84, WGS84, "+proj=tpeqd +lat_1=0 +lon_1=0 +lat_2=60 +lon_2=60 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 54032, aeqd, WGS84, WGS84, "+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 65061, poly, clrk66, NAD27, "+proj=poly +lat_0=13.47246635277778 +lon_0=-144.7487507055556 +x_0=50000.00000000001 +y_0=50000.00000000001 +ellps=clrk66 +datum=NAD27 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 65161, poly, GRS80, NAD83, "+proj=poly +lat_0=13.47246635277778 +lon_0=-144.7487507055556 +x_0=50000 +y_0=50000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102001, aea, GRS80, NAD83, "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102002, lcc, GRS80, NAD83, "+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102003, aea, GRS80, NAD83, "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102004, lcc, GRS80, NAD83, "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102005, eqdc, GRS80, NAD83, "+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102006, aea, GRS80, NAD83, "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102007, aea, GRS80, NAD83, "+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-157 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102008, aea, GRS80, NAD83, "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102009, lcc, GRS80, NAD83, "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102010, eqdc, GRS80, NAD83, "+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102011, sinu, WGS84, WGS84, "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102012, lcc, WGS84, WGS84, "+proj=lcc +lat_1=30 +lat_2=62 +lat_0=0 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102013, aea, intl, "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102014, lcc, intl, "+proj=lcc +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102015, lcc, aust_SA, "+proj=lcc +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102016, aeqd, WGS84, WGS84, "+proj=aeqd +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102017, laea, WGS84, WGS84, "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102018, stere, WGS84, WGS84, "+proj=stere +lat_0=90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102019, aeqd, WGS84, WGS84, "+proj=aeqd +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102020, laea, WGS84, WGS84, "+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102021, stere, WGS84, WGS84, "+proj=stere +lat_0=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102022, aea, WGS84, WGS84, "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102023, eqdc, WGS84, WGS84, "+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=20 +lat_2=-23 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102024, lcc, WGS84, WGS84, "+proj=lcc +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102025, aea, WGS84, WGS84, "+proj=aea +lat_1=15 +lat_2=65 +lat_0=30 +lon_0=95 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102026, eqdc, WGS84, WGS84, "+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=15 +lat_2=65 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102027, lcc, WGS84, WGS84, "+proj=lcc +lat_1=15 +lat_2=65 +lat_0=30 +lon_0=95 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102028, aea, WGS84, WGS84, "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102029, eqdc, WGS84, WGS84, "+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102030, lcc, WGS84, WGS84, "+proj=lcc +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102031, eqdc, intl, "+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=43 +lat_2=62 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102032, eqdc, aust_SA, "+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=-5 +lat_2=-42 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102033, aea, aust_SA, "+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102065, krovak, bessel, "+proj=krovak +lat_0=49.5 +lon_0=24.83333333333333 +alpha=30.28813975277778 +k=0.9999 +x_0=0 +y_0=0 +ellps=bessel +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102066, krovak, bessel, "+proj=krovak +lat_0=49.5 +lon_0=42.5 +alpha=30.28813975277778 +k=0.9999 +x_0=0 +y_0=0 +ellps=bessel +pm=-17.66666666666667 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102067, krovak, bessel, "+proj=krovak +lat_0=49.5 +lon_0=24.83333333333333 +alpha=30.28813975277778 +k=0.9999 +x_0=0 +y_0=0 +ellps=bessel +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102091, tmerc, intl, "+proj=tmerc +lat_0=0 +lon_0=9 +k=0.9996 +x_0=1500000 +y_0=0 +ellps=intl +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102092, tmerc, intl, "+proj=tmerc +lat_0=0 +lon_0=15 +k=0.9996 +x_0=2520000 +y_0=0 +ellps=intl +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102101, tmerc, 6377492.018, 6356173.508712696, "+proj=tmerc +lat_0=58 +lon_0=6.05625 +k=1 +x_0=0 +y_0=0 +a=6377492.018 +b=6356173.508712696 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102102, tmerc, 6377492.018, 6356173.508712696, "+proj=tmerc +lat_0=58 +lon_0=8.389583333333333 +k=1 +x_0=0 +y_0=0 +a=6377492.018 +b=6356173.508712696 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102103, tmerc, 6377492.018, 6356173.508712696, "+proj=tmerc +lat_0=58 +lon_0=10.72291666666667 +k=1 +x_0=0 +y_0=0 +a=6377492.018 +b=6356173.508712696 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102104, tmerc, 6377492.018, 6356173.508712696, "+proj=tmerc +lat_0=58 +lon_0=13.22291666666667 +k=1 +x_0=0 +y_0=0 +a=6377492.018 +b=6356173.508712696 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102105, tmerc, 6377492.018, 6356173.508712696, "+proj=tmerc +lat_0=58 +lon_0=16.88958333333333 +k=1 +x_0=0 +y_0=0 +a=6377492.018 +b=6356173.508712696 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102106, tmerc, 6377492.018, 6356173.508712696, "+proj=tmerc +lat_0=58 +lon_0=20.88958333333333 +k=1 +x_0=0 +y_0=0 +a=6377492.018 +b=6356173.508712696 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102107, tmerc, 6377492.018, 6356173.508712696, "+proj=tmerc +lat_0=58 +lon_0=24.88958333333333 +k=1 +x_0=0 +y_0=0 +a=6377492.018 +b=6356173.508712696 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102108, tmerc, 6377492.018, 6356173.508712696, "+proj=tmerc +lat_0=58 +lon_0=29.05625 +k=1 +x_0=0 +y_0=0 +a=6377492.018 +b=6356173.508712696 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102110, lcc, GRS80, "+proj=lcc +lat_1=44 +lat_2=49 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102114, utm, clrk66, "+proj=utm +zone=4 +ellps=clrk66 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102115, utm, clrk66, "+proj=utm +zone=5 +ellps=clrk66 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102120, omerc, clrk66, NAD27, "+proj=omerc +lat_0=45.30916666666666 +lonc=-86 +alpha=337.255555555556 +k=0.9996 +x_0=2546731.495961392 +y_0=-4354009.816002033 +ellps=clrk66 +datum=NAD27 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102121, omerc, GRS80, NAD83, "+proj=omerc +lat_0=45.30916666666666 +lonc=-86 +alpha=337.255555555556 +k=0.9996 +x_0=2546731.495961392 +y_0=-4354009.816002033 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102122, omerc, clrk66, NAD27, "+proj=omerc +lat_0=45.30916666666666 +lonc=-86 +alpha=337.255555555556 +k=0.9996 +x_0=2546731.496 +y_0=-4354009.816 +ellps=clrk66 +datum=NAD27 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102123, omerc, GRS80, NAD83, "+proj=omerc +lat_0=45.30916666666666 +lonc=-86 +alpha=337.255555555556 +k=0.9996 +x_0=2546731.496 +y_0=-4354009.816 +ellps=GRS80 +datum=NAD83 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102132, utm, 6377492.018, 6356173.508712696, "+proj=utm +zone=32 +a=6377492.018 +b=6356173.508712696 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102133, utm, 6377492.018, 6356173.508712696, "+proj=utm +zone=33 +a=6377492.018 +b=6356173.508712696 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102134, utm, 6377492.018, 6356173.508712696, "+proj=utm +zone=34 +a=6377492.018 +b=6356173.508712696 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102135, utm, 6377492.018, 6356173.508712696, "+proj=utm +zone=35 +a=6377492.018 +b=6356173.508712696 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102140, tmerc, intl, "+proj=tmerc +lat_0=22.31213333333334 +lon_0=114.1785555555556 +k=1 +x_0=836694.05 +y_0=819069.8 +ellps=intl +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102141, utm, intl, "+proj=utm +zone=49 +ellps=intl +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102142, utm, intl, "+proj=utm +zone=50 +ellps=intl +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102151, utm, bessel, "+proj=utm +zone=51 +ellps=bessel +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102152, utm, bessel, "+proj=utm +zone=52 +ellps=bessel +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102153, utm, bessel, "+proj=utm +zone=53 +ellps=bessel +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102154, utm, bessel, "+proj=utm +zone=54 +ellps=bessel +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102155, utm, bessel, "+proj=utm +zone=55 +ellps=bessel +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102156, utm, bessel, "+proj=utm +zone=56 +ellps=bessel +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102160, tmerc, intl, "+proj=tmerc +lat_0=39.66666666666666 +lon_0=-8.131906111111112 +k=1 +x_0=200180.598 +y_0=299913.01 +ellps=intl +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102161, tmerc, intl, "+proj=tmerc +lat_0=39.66666666666666 +lon_0=-8.131906111111112 +k=1 +x_0=180.598 +y_0=-86.98999999999999 +ellps=intl +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102162, utm, intl, "+proj=utm +zone=26 +ellps=intl +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102163, bonne, bessel, "+proj=bonne +lon_0=-8.131906111111112 +lat_1=39.66666666666666 +x_0=0 +y_0=0 +ellps=bessel +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102164, tmerc, intl, "+proj=tmerc +lat_0=39.66666666666666 +lon_0=-8.131906111111112 +k=1 +x_0=200000 +y_0=300000 +ellps=intl +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102165, tmerc, intl, "+proj=tmerc +lat_0=39.66666666666666 +lon_0=-8.131906111111112 +k=1 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102166, utm, intl, "+proj=utm +zone=25 +ellps=intl +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102167, utm, intl, "+proj=utm +zone=28 +ellps=intl +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102168, utm, intl, "+proj=utm +zone=26 +ellps=intl +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102169, utm, intl, "+proj=utm +zone=28 +ellps=intl +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102191, lcc, 6378249.2, 6356514.999904194, "+proj=lcc +lat_1=33.3 +lat_0=33.3 +lon_0=-5.4 +k_0=0.999625769 +x_0=500000 +y_0=300000 +a=6378249.2 +b=6356514.999904194 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102192, lcc, 6378249.2, 6356514.999904194, "+proj=lcc +lat_1=29.7 +lat_0=29.7 +lon_0=-5.4 +k_0=0.9996155960000001 +x_0=500000 +y_0=300000 +a=6378249.2 +b=6356514.999904194 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102193, lcc, 6378249.2, 6356514.999904194, "+proj=lcc +lat_1=26.1 +lat_0=26.1 +lon_0=-5.4 +k_0=0.9996 +x_0=1200000 +y_0=400000 +a=6378249.2 +b=6356514.999904194 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102229, tmerc, GRS80, "+proj=tmerc +lat_0=30.5 +lon_0=-85.83333333333333 +k=0.99996 +x_0=200000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102230, tmerc, GRS80, "+proj=tmerc +lat_0=30 +lon_0=-87.5 +k=0.9999333333333333 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102241, lcc, GRS80, "+proj=lcc +lat_1=40 +lat_2=41.66666666666666 +lat_0=39.33333333333334 +lon_0=-122 +x_0=2000000 +y_0=500000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102242, lcc, GRS80, "+proj=lcc +lat_1=38.33333333333334 +lat_2=39.83333333333334 +lat_0=37.66666666666666 +lon_0=-122 +x_0=2000000 +y_0=500000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102243, lcc, GRS80, "+proj=lcc +lat_1=37.06666666666667 +lat_2=38.43333333333333 +lat_0=36.5 +lon_0=-120.5 +x_0=2000000 +y_0=500000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102244, lcc, GRS80, "+proj=lcc +lat_1=36 +lat_2=37.25 +lat_0=35.33333333333334 +lon_0=-119 +x_0=2000000 +y_0=500000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102245, lcc, GRS80, "+proj=lcc +lat_1=34.03333333333333 +lat_2=35.46666666666667 +lat_0=33.5 +lon_0=-118 +x_0=2000000 +y_0=500000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102246, lcc, GRS80, "+proj=lcc +lat_1=32.78333333333333 +lat_2=33.88333333333333 +lat_0=32.16666666666666 +lon_0=-116.25 +x_0=2000000 +y_0=500000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102248, tmerc, GRS80, "+proj=tmerc +lat_0=31 +lon_0=-110.1666666666667 +k=0.9999 +x_0=213360 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102249, tmerc, GRS80, "+proj=tmerc +lat_0=31 +lon_0=-111.9166666666667 +k=0.9999 +x_0=213360 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102250, tmerc, GRS80, "+proj=tmerc +lat_0=31 +lon_0=-113.75 +k=0.9999333333333333 +x_0=213360 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102251, lcc, GRS80, "+proj=lcc +lat_1=34.93333333333333 +lat_2=36.23333333333333 +lat_0=34.33333333333334 +lon_0=-92 +x_0=400000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102252, lcc, GRS80, "+proj=lcc +lat_1=33.3 +lat_2=34.76666666666667 +lat_0=32.66666666666666 +lon_0=-92 +x_0=400000 +y_0=400000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102253, lcc, GRS80, "+proj=lcc +lat_1=39.71666666666667 +lat_2=40.78333333333333 +lat_0=39.33333333333334 +lon_0=-105.5 +x_0=914401.8289 +y_0=304800.6096 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102254, lcc, GRS80, "+proj=lcc +lat_1=38.45 +lat_2=39.75 +lat_0=37.83333333333334 +lon_0=-105.5 +x_0=914401.8289 +y_0=304800.6096 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102255, lcc, GRS80, "+proj=lcc +lat_1=37.23333333333333 +lat_2=38.43333333333333 +lat_0=36.66666666666666 +lon_0=-105.5 +x_0=914401.8289 +y_0=304800.6096 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102256, lcc, GRS80, "+proj=lcc +lat_1=41.2 +lat_2=41.86666666666667 +lat_0=40.83333333333334 +lon_0=-72.75 +x_0=304800.6096 +y_0=152400.3048 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102257, tmerc, GRS80, "+proj=tmerc +lat_0=38 +lon_0=-75.41666666666667 +k=0.999995 +x_0=200000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102258, tmerc, GRS80, "+proj=tmerc +lat_0=24.33333333333333 +lon_0=-81 +k=0.9999411764705882 +x_0=200000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102259, tmerc, GRS80, "+proj=tmerc +lat_0=24.33333333333333 +lon_0=-82 +k=0.9999411764705882 +x_0=200000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102260, lcc, GRS80, "+proj=lcc +lat_1=29.58333333333333 +lat_2=30.75 +lat_0=29 +lon_0=-84.5 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102261, tmerc, GRS80, "+proj=tmerc +lat_0=18.83333333333333 +lon_0=-155.5 +k=0.9999666666666667 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102262, tmerc, GRS80, "+proj=tmerc +lat_0=20.33333333333333 +lon_0=-156.6666666666667 +k=0.9999666666666667 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102263, tmerc, GRS80, "+proj=tmerc +lat_0=21.16666666666667 +lon_0=-158 +k=0.99999 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102264, tmerc, GRS80, "+proj=tmerc +lat_0=21.83333333333333 +lon_0=-159.5 +k=0.99999 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102265, tmerc, GRS80, "+proj=tmerc +lat_0=21.66666666666667 +lon_0=-160.1666666666667 +k=1 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102266, tmerc, GRS80, "+proj=tmerc +lat_0=30 +lon_0=-82.16666666666667 +k=0.9999 +x_0=200000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102267, tmerc, GRS80, "+proj=tmerc +lat_0=30 +lon_0=-84.16666666666667 +k=0.9999 +x_0=700000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102268, tmerc, GRS80, "+proj=tmerc +lat_0=41.66666666666666 +lon_0=-112.1666666666667 +k=0.9999473684210526 +x_0=200000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102269, tmerc, GRS80, "+proj=tmerc +lat_0=41.66666666666666 +lon_0=-114 +k=0.9999473684210526 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102270, tmerc, GRS80, "+proj=tmerc +lat_0=41.66666666666666 +lon_0=-115.75 +k=0.9999333333333333 +x_0=800000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102271, tmerc, GRS80, "+proj=tmerc +lat_0=36.66666666666666 +lon_0=-88.33333333333333 +k=0.9999749999999999 +x_0=300000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102272, tmerc, GRS80, "+proj=tmerc +lat_0=36.66666666666666 +lon_0=-90.16666666666667 +k=0.9999411764705882 +x_0=700000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102273, tmerc, GRS80, "+proj=tmerc +lat_0=37.5 +lon_0=-85.66666666666667 +k=0.9999666666666667 +x_0=100000 +y_0=250000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102274, tmerc, GRS80, "+proj=tmerc +lat_0=37.5 +lon_0=-87.08333333333333 +k=0.9999666666666667 +x_0=900000 +y_0=250000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102277, lcc, GRS80, "+proj=lcc +lat_1=38.71666666666667 +lat_2=39.78333333333333 +lat_0=38.33333333333334 +lon_0=-98 +x_0=400000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102278, lcc, GRS80, "+proj=lcc +lat_1=37.26666666666667 +lat_2=38.56666666666667 +lat_0=36.66666666666666 +lon_0=-98.5 +x_0=400000 +y_0=400000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102279, lcc, GRS80, "+proj=lcc +lat_1=37.96666666666667 +lat_2=38.96666666666667 +lat_0=37.5 +lon_0=-84.25 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102280, lcc, GRS80, "+proj=lcc +lat_1=36.73333333333333 +lat_2=37.93333333333333 +lat_0=36.33333333333334 +lon_0=-85.75 +x_0=500000 +y_0=500000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102281, lcc, GRS80, "+proj=lcc +lat_1=31.16666666666667 +lat_2=32.66666666666666 +lat_0=30.5 +lon_0=-92.5 +x_0=1000000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102282, lcc, GRS80, "+proj=lcc +lat_1=29.3 +lat_2=30.7 +lat_0=28.5 +lon_0=-91.33333333333333 +x_0=1000000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102283, tmerc, GRS80, "+proj=tmerc +lat_0=43.66666666666666 +lon_0=-68.5 +k=0.9999 +x_0=300000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102284, tmerc, GRS80, "+proj=tmerc +lat_0=42.83333333333334 +lon_0=-70.16666666666667 +k=0.9999666666666667 +x_0=900000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102285, lcc, GRS80, "+proj=lcc +lat_1=38.3 +lat_2=39.45 +lat_0=37.66666666666666 +lon_0=-77 +x_0=400000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102286, lcc, GRS80, "+proj=lcc +lat_1=41.71666666666667 +lat_2=42.68333333333333 +lat_0=41 +lon_0=-71.5 +x_0=200000 +y_0=750000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102287, lcc, GRS80, "+proj=lcc +lat_1=41.28333333333333 +lat_2=41.48333333333333 +lat_0=41 +lon_0=-70.5 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102288, lcc, GRS80, "+proj=lcc +lat_1=45.48333333333333 +lat_2=47.08333333333334 +lat_0=44.78333333333333 +lon_0=-87 +x_0=8000000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102289, lcc, GRS80, "+proj=lcc +lat_1=44.18333333333333 +lat_2=45.7 +lat_0=43.31666666666667 +lon_0=-84.36666666666666 +x_0=6000000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102290, lcc, GRS80, "+proj=lcc +lat_1=42.1 +lat_2=43.66666666666666 +lat_0=41.5 +lon_0=-84.36666666666666 +x_0=4000000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102291, lcc, GRS80, "+proj=lcc +lat_1=47.03333333333333 +lat_2=48.63333333333333 +lat_0=46.5 +lon_0=-93.09999999999999 +x_0=800000 +y_0=100000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102292, lcc, GRS80, "+proj=lcc +lat_1=45.61666666666667 +lat_2=47.05 +lat_0=45 +lon_0=-94.25 +x_0=800000 +y_0=100000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102293, lcc, GRS80, "+proj=lcc +lat_1=43.78333333333333 +lat_2=45.21666666666667 +lat_0=43 +lon_0=-94 +x_0=800000 +y_0=100000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102294, tmerc, GRS80, "+proj=tmerc +lat_0=29.5 +lon_0=-88.83333333333333 +k=0.99995 +x_0=300000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102295, tmerc, GRS80, "+proj=tmerc +lat_0=29.5 +lon_0=-90.33333333333333 +k=0.99995 +x_0=700000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102296, tmerc, GRS80, "+proj=tmerc +lat_0=35.83333333333334 +lon_0=-90.5 +k=0.9999333333333333 +x_0=250000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102297, tmerc, GRS80, "+proj=tmerc +lat_0=35.83333333333334 +lon_0=-92.5 +k=0.9999333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102298, tmerc, GRS80, "+proj=tmerc +lat_0=36.16666666666666 +lon_0=-94.5 +k=0.9999411764705882 +x_0=850000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102300, lcc, GRS80, "+proj=lcc +lat_1=45 +lat_2=49 +lat_0=44.25 +lon_0=-109.5 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102304, lcc, GRS80, "+proj=lcc +lat_1=40 +lat_2=43 +lat_0=39.83333333333334 +lon_0=-100 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102307, tmerc, GRS80, "+proj=tmerc +lat_0=34.75 +lon_0=-115.5833333333333 +k=0.9999 +x_0=200000 +y_0=8000000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102308, tmerc, GRS80, "+proj=tmerc +lat_0=34.75 +lon_0=-116.6666666666667 +k=0.9999 +x_0=500000 +y_0=6000000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102309, tmerc, GRS80, "+proj=tmerc +lat_0=34.75 +lon_0=-118.5833333333333 +k=0.9999 +x_0=800000 +y_0=4000000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102310, tmerc, GRS80, "+proj=tmerc +lat_0=42.5 +lon_0=-71.66666666666667 +k=0.9999666666666667 +x_0=300000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102311, tmerc, GRS80, "+proj=tmerc +lat_0=38.83333333333334 +lon_0=-74.5 +k=0.9999 +x_0=150000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102312, tmerc, GRS80, "+proj=tmerc +lat_0=31 +lon_0=-104.3333333333333 +k=0.9999090909090909 +x_0=165000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102313, tmerc, GRS80, "+proj=tmerc +lat_0=31 +lon_0=-106.25 +k=0.9999 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102314, tmerc, GRS80, "+proj=tmerc +lat_0=31 +lon_0=-107.8333333333333 +k=0.9999166666666667 +x_0=830000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102315, tmerc, GRS80, "+proj=tmerc +lat_0=38.83333333333334 +lon_0=-74.5 +k=0.9999 +x_0=150000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102316, tmerc, GRS80, "+proj=tmerc +lat_0=40 +lon_0=-76.58333333333333 +k=0.9999375 +x_0=250000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102317, tmerc, GRS80, "+proj=tmerc +lat_0=40 +lon_0=-78.58333333333333 +k=0.9999375 +x_0=350000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102318, lcc, GRS80, "+proj=lcc +lat_1=40.66666666666666 +lat_2=41.03333333333333 +lat_0=40.16666666666666 +lon_0=-74 +x_0=300000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102320, lcc, GRS80, "+proj=lcc +lat_1=47.43333333333333 +lat_2=48.73333333333333 +lat_0=47 +lon_0=-100.5 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102321, lcc, GRS80, "+proj=lcc +lat_1=46.18333333333333 +lat_2=47.48333333333333 +lat_0=45.66666666666666 +lon_0=-100.5 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102322, lcc, GRS80, "+proj=lcc +lat_1=40.43333333333333 +lat_2=41.7 +lat_0=39.66666666666666 +lon_0=-82.5 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102323, lcc, GRS80, "+proj=lcc +lat_1=38.73333333333333 +lat_2=40.03333333333333 +lat_0=38 +lon_0=-82.5 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102324, lcc, GRS80, "+proj=lcc +lat_1=35.56666666666667 +lat_2=36.76666666666667 +lat_0=35 +lon_0=-98 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102325, lcc, GRS80, "+proj=lcc +lat_1=33.93333333333333 +lat_2=35.23333333333333 +lat_0=33.33333333333334 +lon_0=-98 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102326, lcc, GRS80, "+proj=lcc +lat_1=44.33333333333334 +lat_2=46 +lat_0=43.66666666666666 +lon_0=-120.5 +x_0=2500000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102327, lcc, GRS80, "+proj=lcc +lat_1=42.33333333333334 +lat_2=44 +lat_0=41.66666666666666 +lon_0=-120.5 +x_0=1500000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102330, tmerc, GRS80, "+proj=tmerc +lat_0=41.08333333333334 +lon_0=-71.5 +k=0.99999375 +x_0=100000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102334, lcc, GRS80, "+proj=lcc +lat_1=44.41666666666666 +lat_2=45.68333333333333 +lat_0=43.83333333333334 +lon_0=-100 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102335, lcc, GRS80, "+proj=lcc +lat_1=42.83333333333334 +lat_2=44.4 +lat_0=42.33333333333334 +lon_0=-100.3333333333333 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102336, lcc, GRS80, "+proj=lcc +lat_1=35.25 +lat_2=36.41666666666666 +lat_0=34.33333333333334 +lon_0=-86 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102337, lcc, GRS80, "+proj=lcc +lat_1=34.65 +lat_2=36.18333333333333 +lat_0=34 +lon_0=-101.5 +x_0=200000 +y_0=1000000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102338, lcc, GRS80, "+proj=lcc +lat_1=32.13333333333333 +lat_2=33.96666666666667 +lat_0=31.66666666666667 +lon_0=-98.5 +x_0=600000 +y_0=2000000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102339, lcc, GRS80, "+proj=lcc +lat_1=30.11666666666667 +lat_2=31.88333333333333 +lat_0=29.66666666666667 +lon_0=-100.3333333333333 +x_0=700000 +y_0=3000000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102340, lcc, GRS80, "+proj=lcc +lat_1=28.38333333333333 +lat_2=30.28333333333334 +lat_0=27.83333333333333 +lon_0=-99 +x_0=600000 +y_0=4000000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102341, lcc, GRS80, "+proj=lcc +lat_1=26.16666666666667 +lat_2=27.83333333333333 +lat_0=25.66666666666667 +lon_0=-98.5 +x_0=300000 +y_0=5000000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102342, lcc, GRS80, "+proj=lcc +lat_1=40.71666666666667 +lat_2=41.78333333333333 +lat_0=40.33333333333334 +lon_0=-111.5 +x_0=500000 +y_0=1000000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102343, lcc, GRS80, "+proj=lcc +lat_1=39.01666666666667 +lat_2=40.65 +lat_0=38.33333333333334 +lon_0=-111.5 +x_0=500000 +y_0=2000000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102344, lcc, GRS80, "+proj=lcc +lat_1=37.21666666666667 +lat_2=38.35 +lat_0=36.66666666666666 +lon_0=-111.5 +x_0=500000 +y_0=3000000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102345, tmerc, GRS80, "+proj=tmerc +lat_0=42.5 +lon_0=-72.5 +k=0.9999642857142857 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102346, lcc, GRS80, "+proj=lcc +lat_1=38.03333333333333 +lat_2=39.2 +lat_0=37.66666666666666 +lon_0=-78.5 +x_0=3500000 +y_0=2000000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102347, lcc, GRS80, "+proj=lcc +lat_1=36.76666666666667 +lat_2=37.96666666666667 +lat_0=36.33333333333334 +lon_0=-78.5 +x_0=3500000 +y_0=1000000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102348, lcc, GRS80, "+proj=lcc +lat_1=47.5 +lat_2=48.73333333333333 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102349, lcc, GRS80, "+proj=lcc +lat_1=45.83333333333334 +lat_2=47.33333333333334 +lat_0=45.33333333333334 +lon_0=-120.5 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102350, lcc, GRS80, "+proj=lcc +lat_1=39 +lat_2=40.25 +lat_0=38.5 +lon_0=-79.5 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102351, lcc, GRS80, "+proj=lcc +lat_1=37.48333333333333 +lat_2=38.88333333333333 +lat_0=37 +lon_0=-81 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102352, lcc, GRS80, "+proj=lcc +lat_1=45.56666666666667 +lat_2=46.76666666666667 +lat_0=45.16666666666666 +lon_0=-90 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102353, lcc, GRS80, "+proj=lcc +lat_1=44.25 +lat_2=45.5 +lat_0=43.83333333333334 +lon_0=-90 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102354, lcc, GRS80, "+proj=lcc +lat_1=42.73333333333333 +lat_2=44.06666666666667 +lat_0=42 +lon_0=-90 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102355, tmerc, GRS80, "+proj=tmerc +lat_0=40.5 +lon_0=-105.1666666666667 +k=0.9999375 +x_0=200000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102356, tmerc, GRS80, "+proj=tmerc +lat_0=40.5 +lon_0=-107.3333333333333 +k=0.9999375 +x_0=400000 +y_0=100000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102357, tmerc, GRS80, "+proj=tmerc +lat_0=40.5 +lon_0=-108.75 +k=0.9999375 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102358, tmerc, GRS80, "+proj=tmerc +lat_0=40.5 +lon_0=-110.0833333333333 +k=0.9999375 +x_0=800000 +y_0=100000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102361, lcc, GRS80, "+proj=lcc +lat_1=18.03333333333334 +lat_2=18.43333333333333 +lat_0=17.83333333333333 +lon_0=-66.43333333333334 +x_0=200000 +y_0=200000 +ellps=GRS80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102491, lcc, 6378249.2, 6356514.999904194, "+proj=lcc +lat_1=36 +lat_0=36 +lon_0=2.7 +k_0=0.999625544 +x_0=500000 +y_0=300000 +a=6378249.2 +b=6356514.999904194 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102492, lcc, 6378249.2, 6356514.999904194, "+proj=lcc +lat_1=33.3 +lat_0=33.3 +lon_0=2.7 +k_0=0.999625769 +x_0=500000 +y_0=300000 +a=6378249.2 +b=6356514.999904194 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102581, lcc, 6378249.2, 6356514.999904194, "+proj=lcc +lat_1=49.5 +lat_0=49.5 +lon_0=2.337229166666667 +k_0=0.999877341 +x_0=600000 +y_0=1200000 +a=6378249.2 +b=6356514.999904194 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102582, lcc, 6378249.2, 6356514.999904194, "+proj=lcc +lat_1=46.8 +lat_0=46.8 +lon_0=2.337229166666667 +k_0=0.99987742 +x_0=600000 +y_0=2200000 +a=6378249.2 +b=6356514.999904194 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102583, lcc, 6378249.2, 6356514.999904194, "+proj=lcc +lat_1=44.1 +lat_0=44.1 +lon_0=2.337229166666667 +k_0=0.999877499 +x_0=600000 +y_0=3200000 +a=6378249.2 +b=6356514.999904194 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 102584, lcc, 6378249.2, 6356514.999904194, "+proj=lcc +lat_1=42.165 +lat_0=42.165 +lon_0=2.337229166666667 +k_0=0.99994471 +x_0=234.358 +y_0=4185861.369 +a=6378249.2 +b=6356514.999904194 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102591, lcc, clrk80, "+proj=lcc +lat_1=36 +lat_0=36 +lon_0=2.7 +k_0=0.999625544 +x_0=500135 +y_0=300090 +ellps=clrk80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 102592, lcc, clrk80, "+proj=lcc +lat_1=33.3 +lat_0=33.3 +lon_0=2.7 +k_0=0.999625769 +x_0=500135 +y_0=300090 +ellps=clrk80 +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102629, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=30.5 +lon_0=-85.83333333333333 +k=0.99996 +x_0=200000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102630, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=30 +lon_0=-87.5 +k=0.9999333333333333 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102631, omerc, GRS80, NAD83, "+proj=omerc +lat_0=57 +lonc=-133.6666666666667 +alpha=-36.86989764583333 +k=0.9999 +x_0=4999999.999999999 +y_0=-4999999.999999999 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102632, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=54 +lon_0=-142 +k=0.9999 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102633, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=54 +lon_0=-146 +k=0.9999 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102634, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=54 +lon_0=-150 +k=0.9999 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102635, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=54 +lon_0=-154 +k=0.9999 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102636, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=54 +lon_0=-158 +k=0.9999 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102637, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=54 +lon_0=-162 +k=0.9999 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102638, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=54 +lon_0=-166 +k=0.9999 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102639, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=54 +lon_0=-170 +k=0.9999 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102640, lcc, GRS80, NAD83, "+proj=lcc +lat_1=51.83333333333334 +lat_2=53.83333333333334 +lat_0=51 +lon_0=-176 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102641, lcc, GRS80, NAD83, "+proj=lcc +lat_1=40 +lat_2=41.66666666666666 +lat_0=39.33333333333334 +lon_0=-122 +x_0=2000000 +y_0=500000.0000000002 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102642, lcc, GRS80, NAD83, "+proj=lcc +lat_1=38.33333333333334 +lat_2=39.83333333333334 +lat_0=37.66666666666666 +lon_0=-122 +x_0=2000000 +y_0=500000.0000000002 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102643, lcc, GRS80, NAD83, "+proj=lcc +lat_1=37.06666666666667 +lat_2=38.43333333333333 +lat_0=36.5 +lon_0=-120.5 +x_0=2000000 +y_0=500000.0000000002 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102644, lcc, GRS80, NAD83, "+proj=lcc +lat_1=36 +lat_2=37.25 +lat_0=35.33333333333334 +lon_0=-119 +x_0=2000000 +y_0=500000.0000000002 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102645, lcc, GRS80, NAD83, "+proj=lcc +lat_1=34.03333333333333 +lat_2=35.46666666666667 +lat_0=33.5 +lon_0=-118 +x_0=2000000 +y_0=500000.0000000002 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102646, lcc, GRS80, NAD83, "+proj=lcc +lat_1=32.78333333333333 +lat_2=33.88333333333333 +lat_0=32.16666666666666 +lon_0=-116.25 +x_0=2000000 +y_0=500000.0000000002 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102648, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=31 +lon_0=-110.1666666666667 +k=0.9999 +x_0=213360 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102649, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=31 +lon_0=-111.9166666666667 +k=0.9999 +x_0=213360 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102650, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=31 +lon_0=-113.75 +k=0.9999333333333333 +x_0=213360 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102651, lcc, GRS80, NAD83, "+proj=lcc +lat_1=34.93333333333333 +lat_2=36.23333333333333 +lat_0=34.33333333333334 +lon_0=-92 +x_0=399999.9999999999 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102652, lcc, GRS80, NAD83, "+proj=lcc +lat_1=33.3 +lat_2=34.76666666666667 +lat_0=32.66666666666666 +lon_0=-92 +x_0=399999.9999999999 +y_0=399999.9999999999 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102653, lcc, GRS80, NAD83, "+proj=lcc +lat_1=39.71666666666667 +lat_2=40.78333333333333 +lat_0=39.33333333333334 +lon_0=-105.5 +x_0=914401.8289 +y_0=304800.6096 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102654, lcc, GRS80, NAD83, "+proj=lcc +lat_1=38.45 +lat_2=39.75 +lat_0=37.83333333333334 +lon_0=-105.5 +x_0=914401.8289 +y_0=304800.6096 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102655, lcc, GRS80, NAD83, "+proj=lcc +lat_1=37.23333333333333 +lat_2=38.43333333333333 +lat_0=36.66666666666666 +lon_0=-105.5 +x_0=914401.8289 +y_0=304800.6096 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102656, lcc, GRS80, NAD83, "+proj=lcc +lat_1=41.2 +lat_2=41.86666666666667 +lat_0=40.83333333333334 +lon_0=-72.75 +x_0=304800.6096 +y_0=152400.3048 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102657, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=38 +lon_0=-75.41666666666667 +k=0.999995 +x_0=200000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102658, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=24.33333333333333 +lon_0=-81 +k=0.9999411764705882 +x_0=200000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102659, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=24.33333333333333 +lon_0=-82 +k=0.9999411764705882 +x_0=200000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102660, lcc, GRS80, NAD83, "+proj=lcc +lat_1=29.58333333333333 +lat_2=30.75 +lat_0=29 +lon_0=-84.5 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102661, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=18.83333333333333 +lon_0=-155.5 +k=0.9999666666666667 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102662, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=20.33333333333333 +lon_0=-156.6666666666667 +k=0.9999666666666667 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102663, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=21.16666666666667 +lon_0=-158 +k=0.99999 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102664, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=21.83333333333333 +lon_0=-159.5 +k=0.99999 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102665, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=21.66666666666667 +lon_0=-160.1666666666667 +k=1 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102666, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=30 +lon_0=-82.16666666666667 +k=0.9999 +x_0=200000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102667, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=30 +lon_0=-84.16666666666667 +k=0.9999 +x_0=700000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102668, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=41.66666666666666 +lon_0=-112.1666666666667 +k=0.9999473684210526 +x_0=200000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102669, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=41.66666666666666 +lon_0=-114 +k=0.9999473684210526 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102670, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=41.66666666666666 +lon_0=-115.75 +k=0.9999333333333333 +x_0=799999.9999999999 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102671, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=36.66666666666666 +lon_0=-88.33333333333333 +k=0.9999749999999999 +x_0=300000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102672, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=36.66666666666666 +lon_0=-90.16666666666667 +k=0.9999411764705882 +x_0=700000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102673, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=37.5 +lon_0=-85.66666666666667 +k=0.9999666666666667 +x_0=100000 +y_0=250000 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102674, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=37.5 +lon_0=-87.08333333333333 +k=0.9999666666666667 +x_0=900000.0000000001 +y_0=250000 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102675, lcc, GRS80, NAD83, "+proj=lcc +lat_1=42.06666666666667 +lat_2=43.26666666666667 +lat_0=41.5 +lon_0=-93.5 +x_0=1500000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102676, lcc, GRS80, NAD83, "+proj=lcc +lat_1=40.61666666666667 +lat_2=41.78333333333333 +lat_0=40 +lon_0=-93.5 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102677, lcc, GRS80, NAD83, "+proj=lcc +lat_1=38.71666666666667 +lat_2=39.78333333333333 +lat_0=38.33333333333334 +lon_0=-98 +x_0=399999.9999999999 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102678, lcc, GRS80, NAD83, "+proj=lcc +lat_1=37.26666666666667 +lat_2=38.56666666666667 +lat_0=36.66666666666666 +lon_0=-98.5 +x_0=399999.9999999999 +y_0=399999.9999999999 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102679, lcc, GRS80, NAD83, "+proj=lcc +lat_1=37.96666666666667 +lat_2=38.96666666666667 +lat_0=37.5 +lon_0=-84.25 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102680, lcc, GRS80, NAD83, "+proj=lcc +lat_1=36.73333333333333 +lat_2=37.93333333333333 +lat_0=36.33333333333334 +lon_0=-85.75 +x_0=500000.0000000002 +y_0=500000.0000000002 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102681, lcc, GRS80, NAD83, "+proj=lcc +lat_1=31.16666666666667 +lat_2=32.66666666666666 +lat_0=30.5 +lon_0=-92.5 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102682, lcc, GRS80, NAD83, "+proj=lcc +lat_1=29.3 +lat_2=30.7 +lat_0=28.5 +lon_0=-91.33333333333333 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102683, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=43.66666666666666 +lon_0=-68.5 +k=0.9999 +x_0=300000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102684, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=42.83333333333334 +lon_0=-70.16666666666667 +k=0.9999666666666667 +x_0=900000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102685, lcc, GRS80, NAD83, "+proj=lcc +lat_1=38.3 +lat_2=39.45 +lat_0=37.66666666666666 +lon_0=-77 +x_0=399999.9999999999 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102686, lcc, GRS80, NAD83, "+proj=lcc +lat_1=41.71666666666667 +lat_2=42.68333333333333 +lat_0=41 +lon_0=-71.5 +x_0=200000 +y_0=750000.0000000001 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102687, lcc, GRS80, NAD83, "+proj=lcc +lat_1=41.28333333333333 +lat_2=41.48333333333333 +lat_0=41 +lon_0=-70.5 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102688, lcc, GRS80, NAD83, "+proj=lcc +lat_1=45.48333333333333 +lat_2=47.08333333333334 +lat_0=44.78333333333333 +lon_0=-87 +x_0=7999999.999999999 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102689, lcc, GRS80, NAD83, "+proj=lcc +lat_1=44.18333333333333 +lat_2=45.7 +lat_0=43.31666666666667 +lon_0=-84.36666666666666 +x_0=6000000.000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102690, lcc, GRS80, NAD83, "+proj=lcc +lat_1=42.1 +lat_2=43.66666666666666 +lat_0=41.5 +lon_0=-84.36666666666666 +x_0=4000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102691, lcc, GRS80, NAD83, "+proj=lcc +lat_1=47.03333333333333 +lat_2=48.63333333333333 +lat_0=46.5 +lon_0=-93.09999999999999 +x_0=799999.9999999999 +y_0=100000 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102692, lcc, GRS80, NAD83, "+proj=lcc +lat_1=45.61666666666667 +lat_2=47.05 +lat_0=45 +lon_0=-94.25 +x_0=799999.9999999999 +y_0=100000 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102693, lcc, GRS80, NAD83, "+proj=lcc +lat_1=43.78333333333333 +lat_2=45.21666666666667 +lat_0=43 +lon_0=-94 +x_0=799999.9999999999 +y_0=100000 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102694, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=29.5 +lon_0=-88.83333333333333 +k=0.99995 +x_0=300000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102695, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=29.5 +lon_0=-90.33333333333333 +k=0.99995 +x_0=700000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102696, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=35.83333333333334 +lon_0=-90.5 +k=0.9999333333333333 +x_0=250000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102697, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=35.83333333333334 +lon_0=-92.5 +k=0.9999333333333333 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102698, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=36.16666666666666 +lon_0=-94.5 +k=0.9999411764705882 +x_0=850000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102700, lcc, GRS80, NAD83, "+proj=lcc +lat_1=45 +lat_2=49 +lat_0=44.25 +lon_0=-109.5 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102704, lcc, GRS80, NAD83, "+proj=lcc +lat_1=40 +lat_2=43 +lat_0=39.83333333333334 +lon_0=-100 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102707, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=34.75 +lon_0=-115.5833333333333 +k=0.9999 +x_0=200000 +y_0=7999999.999999999 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102708, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=34.75 +lon_0=-116.6666666666667 +k=0.9999 +x_0=500000.0000000002 +y_0=6000000.000000001 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102709, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=34.75 +lon_0=-118.5833333333333 +k=0.9999 +x_0=799999.9999999999 +y_0=4000000 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102710, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=42.5 +lon_0=-71.66666666666667 +k=0.9999666666666667 +x_0=300000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102711, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=38.83333333333334 +lon_0=-74.5 +k=0.9999 +x_0=150000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102712, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=31 +lon_0=-104.3333333333333 +k=0.9999090909090909 +x_0=165000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102713, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=31 +lon_0=-106.25 +k=0.9999 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102714, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=31 +lon_0=-107.8333333333333 +k=0.9999166666666667 +x_0=829999.9999999999 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102715, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=38.83333333333334 +lon_0=-74.5 +k=0.9999 +x_0=150000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102716, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=40 +lon_0=-76.58333333333333 +k=0.9999375 +x_0=250000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102717, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=40 +lon_0=-78.58333333333333 +k=0.9999375 +x_0=350000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102718, lcc, GRS80, NAD83, "+proj=lcc +lat_1=40.66666666666666 +lat_2=41.03333333333333 +lat_0=40.16666666666666 +lon_0=-74 +x_0=300000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102719, lcc, GRS80, NAD83, "+proj=lcc +lat_1=34.33333333333334 +lat_2=36.16666666666666 +lat_0=33.75 +lon_0=-79 +x_0=609601.2199999999 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102720, lcc, GRS80, NAD83, "+proj=lcc +lat_1=47.43333333333333 +lat_2=48.73333333333333 +lat_0=47 +lon_0=-100.5 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102721, lcc, GRS80, NAD83, "+proj=lcc +lat_1=46.18333333333333 +lat_2=47.48333333333333 +lat_0=45.66666666666666 +lon_0=-100.5 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102722, lcc, GRS80, NAD83, "+proj=lcc +lat_1=40.43333333333333 +lat_2=41.7 +lat_0=39.66666666666666 +lon_0=-82.5 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102723, lcc, GRS80, NAD83, "+proj=lcc +lat_1=38.73333333333333 +lat_2=40.03333333333333 +lat_0=38 +lon_0=-82.5 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102724, lcc, GRS80, NAD83, "+proj=lcc +lat_1=35.56666666666667 +lat_2=36.76666666666667 +lat_0=35 +lon_0=-98 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102725, lcc, GRS80, NAD83, "+proj=lcc +lat_1=33.93333333333333 +lat_2=35.23333333333333 +lat_0=33.33333333333334 +lon_0=-98 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102726, lcc, GRS80, NAD83, "+proj=lcc +lat_1=44.33333333333334 +lat_2=46 +lat_0=43.66666666666666 +lon_0=-120.5 +x_0=2500000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102727, lcc, GRS80, NAD83, "+proj=lcc +lat_1=42.33333333333334 +lat_2=44 +lat_0=41.66666666666666 +lon_0=-120.5 +x_0=1500000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102728, lcc, GRS80, NAD83, "+proj=lcc +lat_1=40.88333333333333 +lat_2=41.95 +lat_0=40.16666666666666 +lon_0=-77.75 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102729, lcc, GRS80, NAD83, "+proj=lcc +lat_1=39.93333333333333 +lat_2=40.96666666666667 +lat_0=39.33333333333334 +lon_0=-77.75 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102730, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=41.08333333333334 +lon_0=-71.5 +k=0.99999375 +x_0=100000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102733, lcc, GRS80, NAD83, "+proj=lcc +lat_1=32.5 +lat_2=34.83333333333334 +lat_0=31.83333333333333 +lon_0=-81 +x_0=609600.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102734, lcc, GRS80, NAD83, "+proj=lcc +lat_1=44.41666666666666 +lat_2=45.68333333333333 +lat_0=43.83333333333334 +lon_0=-100 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102735, lcc, GRS80, NAD83, "+proj=lcc +lat_1=42.83333333333334 +lat_2=44.4 +lat_0=42.33333333333334 +lon_0=-100.3333333333333 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102736, lcc, GRS80, NAD83, "+proj=lcc +lat_1=35.25 +lat_2=36.41666666666666 +lat_0=34.33333333333334 +lon_0=-86 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102737, lcc, GRS80, NAD83, "+proj=lcc +lat_1=34.65 +lat_2=36.18333333333333 +lat_0=34 +lon_0=-101.5 +x_0=200000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102738, lcc, GRS80, NAD83, "+proj=lcc +lat_1=32.13333333333333 +lat_2=33.96666666666667 +lat_0=31.66666666666667 +lon_0=-98.5 +x_0=600000.0000000001 +y_0=2000000 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102739, lcc, GRS80, NAD83, "+proj=lcc +lat_1=30.11666666666667 +lat_2=31.88333333333333 +lat_0=29.66666666666667 +lon_0=-100.3333333333333 +x_0=700000 +y_0=3000000 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102740, lcc, GRS80, NAD83, "+proj=lcc +lat_1=28.38333333333333 +lat_2=30.28333333333334 +lat_0=27.83333333333333 +lon_0=-99 +x_0=600000.0000000001 +y_0=4000000 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102741, lcc, GRS80, NAD83, "+proj=lcc +lat_1=26.16666666666667 +lat_2=27.83333333333333 +lat_0=25.66666666666667 +lon_0=-98.5 +x_0=300000 +y_0=4999999.999999999 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102742, lcc, GRS80, NAD83, "+proj=lcc +lat_1=40.71666666666667 +lat_2=41.78333333333333 +lat_0=40.33333333333334 +lon_0=-111.5 +x_0=500000.0000000002 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102743, lcc, GRS80, NAD83, "+proj=lcc +lat_1=39.01666666666667 +lat_2=40.65 +lat_0=38.33333333333334 +lon_0=-111.5 +x_0=500000.0000000002 +y_0=2000000 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102744, lcc, GRS80, NAD83, "+proj=lcc +lat_1=37.21666666666667 +lat_2=38.35 +lat_0=36.66666666666666 +lon_0=-111.5 +x_0=500000.0000000002 +y_0=3000000 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102745, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=42.5 +lon_0=-72.5 +k=0.9999642857142857 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102746, lcc, GRS80, NAD83, "+proj=lcc +lat_1=38.03333333333333 +lat_2=39.2 +lat_0=37.66666666666666 +lon_0=-78.5 +x_0=3499999.999999999 +y_0=2000000 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102747, lcc, GRS80, NAD83, "+proj=lcc +lat_1=36.76666666666667 +lat_2=37.96666666666667 +lat_0=36.33333333333334 +lon_0=-78.5 +x_0=3499999.999999999 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102748, lcc, GRS80, NAD83, "+proj=lcc +lat_1=47.5 +lat_2=48.73333333333333 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102749, lcc, GRS80, NAD83, "+proj=lcc +lat_1=45.83333333333334 +lat_2=47.33333333333334 +lat_0=45.33333333333334 +lon_0=-120.5 +x_0=500000.0000000002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102750, lcc, GRS80, NAD83, "+proj=lcc +lat_1=39 +lat_2=40.25 +lat_0=38.5 +lon_0=-79.5 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102751, lcc, GRS80, NAD83, "+proj=lcc +lat_1=37.48333333333333 +lat_2=38.88333333333333 +lat_0=37 +lon_0=-81 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102752, lcc, GRS80, NAD83, "+proj=lcc +lat_1=45.56666666666667 +lat_2=46.76666666666667 +lat_0=45.16666666666666 +lon_0=-90 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102753, lcc, GRS80, NAD83, "+proj=lcc +lat_1=44.25 +lat_2=45.5 +lat_0=43.83333333333334 +lon_0=-90 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102754, lcc, GRS80, NAD83, "+proj=lcc +lat_1=42.73333333333333 +lat_2=44.06666666666667 +lat_0=42 +lon_0=-90 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102755, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=40.5 +lon_0=-105.1666666666667 +k=0.9999375 +x_0=200000 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102756, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=40.5 +lon_0=-107.3333333333333 +k=0.9999375 +x_0=399999.9999999999 +y_0=100000 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102757, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=40.5 +lon_0=-108.75 +k=0.9999375 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102758, tmerc, GRS80, NAD83, "+proj=tmerc +lat_0=40.5 +lon_0=-110.0833333333333 +k=0.9999375 +x_0=799999.9999999999 +y_0=100000 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102761, lcc, GRS80, NAD83, "+proj=lcc +lat_1=18.03333333333334 +lat_2=18.43333333333333 +lat_0=17.83333333333333 +lon_0=-66.43333333333334 +x_0=200000 +y_0=200000 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 102766, poly, GRS80, NAD83, "+proj=poly +lat_0=13.47246635277778 +lon_0=-144.7487507055556 +x_0=49999.99999999999 +y_0=49999.99999999999 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 103300, lcc, intl, "+proj=lcc +lat_1=49.8333339 +lat_2=51.16666733333333 +lat_0=90 +lon_0=4.367486666666666 +x_0=150000.01256 +y_0=5400088.4378 +ellps=intl +units=m +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_ED(esri, 104000, longlat, clrk66, NAD27, "+proj=longlat +ellps=clrk66 +datum=NAD27 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 104101, longlat, bessel, "+proj=longlat +ellps=bessel +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 104102, longlat, bessel, "+proj=longlat +ellps=bessel +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 104103, longlat, clrk80, "+proj=longlat +ellps=clrk80 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 104104, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 104105, longlat, bessel, "+proj=longlat +ellps=bessel +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 104106, longlat, intl, "+proj=longlat +ellps=intl +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 104107, longlat, GRS80, "+proj=longlat +ellps=GRS80 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 104108, longlat, GRS80, "+proj=longlat +ellps=GRS80 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 104261, longlat, 6378249.2, 6356514.999904194, "+proj=longlat +a=6378249.2 +b=6356514.999904194 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_AB(esri, 104304, longlat, 6378249.2, 6356514.999904194, "+proj=longlat +a=6378249.2 +b=6356514.999904194 +no_defs");
BOOST_GEOMETRY_PROJECTIONS_DETAIL_SRID_TRAITS_E (esri, 104305, longlat, clrk80, "+proj=longlat +ellps=clrk80 +no_defs");

} // namespace detail
#endif // DOXYGEN_NO_DETAIL

}}} // namespace boost::geometry::projections


#endif
