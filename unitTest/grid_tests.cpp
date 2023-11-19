#include <boost/test/unit_test.hpp>

#include "../grid.h"
#include "../constants.h"

/**
 * @brief Unit tests for the Grid class
 */
BOOST_AUTO_TEST_SUITE(grid_tests)

BOOST_AUTO_TEST_CASE(grid_test){
    constexpr double leftBound = -10.0;
    constexpr double rightBound = 430.3;
    std::size_t dim = 27;
    Grid grid(leftBound, rightBound, dim);

    BOOST_CHECK_SMALL(grid.getLeftBound() - leftBound, epsilon);
    BOOST_CHECK_SMALL(grid.getRightBound() - rightBound, epsilon);
    BOOST_CHECK(grid.getDim() == dim);
    BOOST_CHECK_SMALL(grid.getNthPoint(0) - leftBound, epsilon);
    BOOST_CHECK_SMALL(grid.getNthPoint(dim - 1) - rightBound, epsilon);

    const double step = (rightBound - leftBound)/(dim - 1.0);
    BOOST_CHECK_SMALL(grid.getNthPoint(dim - 8) - grid.getNthPoint(dim - 10) - 2*step, epsilon);
}
BOOST_AUTO_TEST_SUITE_END()