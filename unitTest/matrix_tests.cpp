#include <boost/test/unit_test.hpp>

#include <cmath>
#include <complex>
#include <vector>

#include "test_utils.h"
#include "../matrix.h"
#include "../hermitianMatrix.h"
#include "../random.h"
#include "../constants.h"

/**
 * @brief  Test basically all matrix functions
 */
BOOST_AUTO_TEST_SUITE(matrix_tests)
// Test basic operations between matrices
BOOST_AUTO_TEST_CASE(matrix_operations_test)
{
    constexpr std::size_t nRows = 100;
    constexpr std::size_t nColumns = 200;
    // Test 1) Sum and subtract two random 100x200 matrices
    Matrix<std::complex<double>> m1 = randomComplexMatrix(nRows, nColumns);
    const Matrix<std::complex<double>> m2 = randomComplexMatrix(nRows, nColumns);
    BOOST_CHECK(m2.getDimX() == m1.getDimX());
    BOOST_CHECK(m2.getDimY() == m1.getDimY());

    Matrix<std::complex<double>> m3 = m1 + m2;
    Matrix<std::complex<double>> m4 = m1 - m2;

    BOOST_CHECK(m3.getDimX() == m1.getDimX());
    BOOST_CHECK(m3.getDimY() == m1.getDimY());
    BOOST_CHECK(m4.getDimX() == m1.getDimX());
    BOOST_CHECK(m4.getDimY() == m1.getDimY());
    for(size_t i = 0; i < m1.getDimX(); i++){
        for(size_t j = 0; j < m1.getDimY(); j++){
            CompareComplexNumbers(m3(i,j), m1(i,j) + m2(i,j));
            CompareComplexNumbers(m4(i,j), m1(i,j) - m2(i,j));
        }
    }

    // Test 2) Check += operator
    m1 += m2;
    CompareMatrices(m1, m3);

    // Test 3) Check - operator;
    Matrix<std::complex<double>> mEmpty =  Matrix<std::complex<double>>(nRows, nColumns);
    CompareMatrices(-m1 + m1, mEmpty);

    // Test 4) Check multiplication by scalar
    const std::complex<double> factor(1.4, 5.22);
    Matrix<std::complex<double>> m5 = m4 * factor;
    BOOST_CHECK(m5.getDimX() == m4.getDimX());
    BOOST_CHECK(m5.getDimY() == m4.getDimY());
    for(size_t i = 0; i < m5.getDimX(); i++){
        for(size_t j = 0; j < m5.getDimY(); j++){
            CompareComplexNumbers(m4(i,j) * factor, m5(i,j));
        }
    }
}

BOOST_AUTO_TEST_CASE(hermitian_matrix_tests){
    //TEST 1) Build the matrix (1, i, -i, 2)
    // and check that the spectrum is the expected one
    HermitianMatrix m1Re(2);
    HermitianMatrix m1Im(2);
    m1Re(0,0) = 1;
    m1Re(0,1) = 0;
    m1Re(1,0) = 0;
    m1Re(1,1) = 2;
    m1Im(0,0) = std::complex<double>(0,0);
    m1Im(0,1) = std::complex<double>(0,1);
    m1Im(1,0) = std::complex<double>(0,-1);
    m1Im(1,1) = std::complex<double>(0,0);

    // sum the two matrices
    m1Re += m1Im;
    HermitianMatrix m1Copy = m1Re;

    // find eigenvalues (this call with destroy the content of m1Re)
    BOOST_CHECK(m1Re.findSpectrumAlg2('N'));
    // eigenvalues must be: 1/2(3 \pm sqrt(5))
    BOOST_CHECK_SMALL(m1Re.getEigenValues()[0] - 0.5*(3-sqrt(5)) , epsilon);
    BOOST_CHECK_SMALL(m1Re.getEigenValues()[1] - 0.5*(3+sqrt(5)) , epsilon);

    // At this point the original m1 is destroyed, we have only eigenvalues but not eigenvectors
    // and we cannot recall findSpectrum
    BOOST_CHECK(!m1Re.findSpectrumAlg2('V'));
    // If we try to find again eigenvalues of course this will just return true, since we already calculated them.
    BOOST_CHECK(m1Re.findSpectrumAlg2('N'));
    // To actually find eigenvectors let's use the copy
    BOOST_CHECK(m1Copy.findSpectrumAlg2('V'));
    // normalize eigenvectors in such a way that the second component is 1
    std::vector<std::complex<double>> eigenVector1 = {m1Copy(0,0), m1Copy(1,0)};
    std::vector<std::complex<double>> eigenVector2 = {m1Copy(0,1), m1Copy(1,1)};
    eigenVector1[0] /= eigenVector1[1];
    eigenVector1[1] /= eigenVector1[1];

    eigenVector2[0] /= eigenVector2[1];
    eigenVector2[1] /= eigenVector2[1];
    // With this normalization eigenvectors should be:
    // v1 = (-i/2 * (1 + sqrt(5)),1)
    // v2 = (i/2 * (-1 + sqrt(5)),1)
    CompareComplexNumbers(eigenVector1[0], {0, -0.5 * (1+sqrt(5))});
    CompareComplexNumbers(eigenVector1[1], 1);

    CompareComplexNumbers(eigenVector2[0], {0, 0.5 * (-1+sqrt(5))});
    CompareComplexNumbers(eigenVector2[1], 1);

    //TEST 2) Check that eigenvalues are in ascending order
    // check that eigenvalues are normalized, for a random big 100x100 matrix
    constexpr double dim = 100;
    HermitianMatrix m2 = randomComplexLTHermitianMatrix(dim);
    m2.findSpectrumAlg2('V');
    for(int i = 0; i < dim-1; i++){
        BOOST_CHECK(m2.getEigenValues()[i+1] >= m2.getEigenValues()[i]);
    }
    double normalization = 0;
    for(int i = 0; i < dim; i++){
        normalization = 0;
        for( int j = 0; j < dim; j++){
            normalization += std::norm(m2.getEigenVectors()(j,i));
        }
        BOOST_CHECK_SMALL(normalization - 1, epsilon);
    }
}

BOOST_AUTO_TEST_CASE(matrix_tensor_product){
    //Test 1:
    // tensor product between
    // (1, 1)  and (1, 0)
    // (1, 1)      (0, 1)
    Matrix<std::complex<double>> m1(2, 2);
    m1(0,0) = 1.0;
    m1(0,1) = 1.0;
    m1(1,0) = 1.0;
    m1(1,1) = 1.0;
    Matrix<std::complex<double>> m2(2, 2);
    m2(0,0) = 1.0;
    m2(1,1) = 1.0;
    // (1, 0, 1, 0)
    // (0, 1, 0, 1)
    // (1, 0, 1, 0)
    // (0, 1, 0, 1)
    Matrix<std::complex<double>> m3 = tens_product(m1,m2);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            CompareComplexNumbers(m3(i, j), ((i + j)%2 == 0) ? 1.0 : 0.0);
        }
    }
    // (1, 1, 0, 0)
    // (1, 1, 0, 0)
    // (0, 0, 1, 1)
    // (0, 0, 1, 1)
    Matrix<std::complex<double>> m4 = tens_product(m2,m1);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            if((i < 2 && j < 2 )|| (i>=2 && j>=2)){
                CompareComplexNumbers(m4(i, j), 1.0);
            }else{
                CompareComplexNumbers(m4(i, j), 0.0);
            }
        }
    }
    // Test 2:
    // tensor product between
    // (2)  (0, 5)
    // (0)
    Matrix<std::complex<double>> m5(2, 1);
    m5(0, 0) = 2.0;
    m5(1, 0) = 0.0;
    Matrix<std::complex<double>> m6(1, 2);
    m6(0,0) = 0.0;
    m6(0,1) = 5.0;
    // (0, 10)
    // (0,  0)
    Matrix<std::complex<double>> m7 = tens_product(m5,m6);
    Matrix<std::complex<double>> m8 = tens_product(m6,m5);
    
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            CompareComplexNumbers(m7(i, j), (i == 0 && j == 1) ? 10.0 : 0.0);
            CompareComplexNumbers(m8(i, j), (i == 0 && j == 1) ? 10.0 : 0.0);
        }
    }
}

BOOST_AUTO_TEST_CASE(matrix_multiplication_tests){
    //Test 1:
    // multiply by identity and verify that result does not change
    // (1, 1)  and (1, 0)
    // (1, 1)      (0, 1)
    Matrix<std::complex<double>> m1(2, 2);
    m1(0,0) = 1.0;
    m1(0,1) = 1.0;
    m1(1,0) = 1.0;
    m1(1,1) = 1.0;
    Matrix<std::complex<double>> m2(2, 2);
    m2(0,0) = 1.0;
    m2(1,1) = 1.0;
    CompareMatrices(m1, matrix_product(m1,m2));
    CompareMatrices(m1, matrix_product(m2,m1));
    // Test 2:
    // matrix product between
    // (2)  (0, 10, 20)
    // (3)
    // (4)
    Matrix<std::complex<double>> m3(3, 1);
    m3(0, 0) = 2.0;
    m3(1, 0) = 3.0;
    m3(2, 0) = 4.0;
    Matrix<std::complex<double>> m4(1, 3);
    m4(0,0) = 0.0;
    m4(0,1) = 10.0;
    m4(0,2) = 20.0;
    // (30 + 80) = 110
    Matrix<std::complex<double>> m5(1, 1);
    m5(0,0) = 110;
    CompareMatrices(m5, matrix_product(m4,m3));
    //Test 3:
    // matrix product between
    // (1 )  (0, 10, 20)
    // (-1)  (30, 40, 50)
    // (2 )
    Matrix<std::complex<double>> m6(3, 1);
    m6(0, 0) = 1.0;
    m6(1, 0) = -1.0;
    m6(2, 0) = 2.0;
    Matrix<std::complex<double>> m7(2, 3);
    m7(0, 0) = 0.0;
    m7(0, 1) = 10.0;
    m7(0, 2) = 20.0;
    m7(1, 0) = 30.0;
    m7(1, 1) = 40.0;
    m7(1, 2) = 50.0;
    Matrix<std::complex<double>> m8(2, 1);
    // (30)
    // (90)
    m8(0, 0) = 30.0;
    m8(1, 0) = 90.0;
    CompareMatrices(m8, matrix_product(m7,m6));
}
BOOST_AUTO_TEST_SUITE_END()