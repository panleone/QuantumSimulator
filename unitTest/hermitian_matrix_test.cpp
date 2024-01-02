#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <complex>
#include <cstddef>
#include <math.h>

#include "test_utils.h"
#include "../hermitianMatrix.h"

/**
 * @brief Unit tests for the Hamiltonian class
 */
BOOST_AUTO_TEST_SUITE(hermitian_matrix_tests)
BOOST_AUTO_TEST_CASE(hermitian_matrix_spectrum){
    constexpr size_t dim = 2;
    // Create the matrix
    // (1 , i)  
    // (-i, 2)
    HermitianMatrix M(2);
    M(0,0) = 1;
    M(0,1) = std::complex<double>(0,1);
    M(1,0) = std::complex<double>(0,-1);
    M(1,1) = 2;
    HermitianMatrix Mcopy{M};
    
    // Find eigenvalues and eigenvectors
    M.findSpectrumAlg2('V');
    for(size_t i = 0; i<2; i++ ){
        //Verify that the ith eigenvector is indeed an eigenvector
        auto ithEigenValue = M.getEigenValues().at(i);
        auto ithEigenVector = matrix_product(Mcopy, M.getNthEigenVector(i));
        ithEigenVector /= ithEigenValue;
        CompareMatrices(M.getNthEigenVector(i), ithEigenVector);
    }
    // Verify that the projector is built correctly
    ComplexVector emptyVector{2};
    for(size_t i = 0; i<2; i++ ){
        auto ithEigenVector = M.getNthEigenVector(i);
        auto otherEigenVector = M.getNthEigenVector((i+1)%2);
        CompareMatrices(emptyVector, matrix_product(ithEigenVector.getProjector(), otherEigenVector ));
        CompareMatrices(ithEigenVector, matrix_product(ithEigenVector.getProjector(), ithEigenVector));
    }
}
BOOST_AUTO_TEST_CASE(adjoint_test){
    constexpr size_t dimX = 2;
    constexpr size_t dimY = 3;
    
    // Create the matrix
    // (1 , i, 2 + 2i)  
    // (0 , 7, 10 + 10i)
    ComplexMatrix M(2,3);
    M(0,0) = 1;
    M(0,1) = std::complex<double>(0,1);
    M(0,2) = std::complex<double>(2,2); 
    M(1,0) = 0.0;
    M(1,1) = 7.0;
    M(1,2) = std::complex<double>(10,10);

    // It's adjoing is the matrix
    // (1     , 0       )  
    // (-i    , 7       )
    // (2 - 2i, 10 - 10i)
    ComplexMatrix adjM(3,2);
    adjM(0,0) = 1;
    adjM(1,0) = std::complex<double>(0,-1);
    adjM(2,0) = std::complex<double>(2,-2); 
    adjM(0,1) = 0.0;
    adjM(1,1) = 7.0;
    adjM(2,1) = std::complex<double>(10,-10);
    CompareMatrices(M.getAdjoint(), adjM);
    CompareMatrices(adjM.getAdjoint(), M);
}
BOOST_AUTO_TEST_SUITE_END()