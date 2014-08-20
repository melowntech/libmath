#include <boost/test/unit_test.hpp>

#include <math/math_all.hpp>

#include "dbglog/dbglog.hpp"

BOOST_AUTO_TEST_CASE(math_extents_scale)
{
    const math::Extents2 e(100, 300, 400, 500);
    const double scale(1.5);

    const auto scaled(e * scale);

    BOOST_REQUIRE(center(scaled) == center(e));
    BOOST_REQUIRE(size(scaled) == (size(e) * scale));
}
