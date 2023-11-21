#ifndef TEST_UTILS
#define TEST_UTILS

#include <boost/test/unit_test.hpp>
#include <complex>
#include <cstddef>

#include "../constants.h"
#include "../matrix.h"

/**
 * @brief BOOST_CHECK that two complex numbers are the same with precision epsilon
 */
inline void CompareComplexNumbers(const std::complex<double>& c1, const std::complex<double>& c2){
    BOOST_CHECK_SMALL(c1.imag() - c2.imag(), epsilon);
    BOOST_CHECK_SMALL(c1.real() - c2.real(), epsilon);
}

/**
 * @brief Check element by element that two matrices are the same, with precision epsilon
 */
inline void CompareMatrices(const Matrix<std::complex<double>>& m1, const Matrix<std::complex<double>>& m2){
    BOOST_CHECK(m1.getDimX() == m2.getDimX());
    BOOST_CHECK(m1.getDimY() == m2.getDimY());
    for(size_t i = 0; i < m1.getDimX(); i++){
        for(size_t j = 0; j < m1.getDimY(); j++){
            CompareComplexNumbers(m1(i,j), m2(i,j));
        }
    }
}

template<typename T>
// Return the sign of T
int sgn(T val){
    return (T(0) < val) - (val < T(0));
}

#endif