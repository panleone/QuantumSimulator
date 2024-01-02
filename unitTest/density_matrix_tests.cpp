#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <cstddef>

#include "../densityMatrix.h"
#include "../multiParticleWF.h"
#include "../constants.h"
#include "test_utils.h"

/**
 * @brief Unit tests for the Grid class
 */
BOOST_AUTO_TEST_SUITE(density_matrix_tests)

// Test that density matrices are constructed in the right way
BOOST_AUTO_TEST_CASE(density_matrix_constructor_test){
    // In the whole test wavefunctions will represent a system of two qubits
    MultiParticleWF waveFunctionTestA(2,2);
    // set it in the state (2+2i)|0>
    waveFunctionTestA(0) = std::complex<double>(2,2);
    waveFunctionTestA(1) = 0;
    waveFunctionTestA(2) = 0;
    waveFunctionTestA(3) = 0;
    BOOST_CHECK(waveFunctionTestA.getDim() == 4);
    BOOST_CHECK_SMALL(waveFunctionTestA.getNorm() - sqrt(8), epsilon);

    // Build the corresponding DensityMatrix:
    // (1, 0, 0, 0)
    // (0, 0, 0, 0)
    // (0, 0, 0, 0)
    // (0, 0, 0, 0)
    PureDensityMatrix densityMatrixA(waveFunctionTestA);
    BOOST_CHECK(densityMatrixA.getDim() == 4);
    BOOST_CHECK(densityMatrixA.getN() == 2);
    CompareComplexNumbers(densityMatrixA.trace(), 1.0);
    for(size_t i = 0; i < 4; i++){
        for(size_t j = 0; j < 4; j++){
            if(i == 0 && j == 0){
                CompareComplexNumbers(densityMatrixA(i,j), 1.0);
            } else {
                CompareComplexNumbers(densityMatrixA(i,j), 0.0);
            }
        }
    }

    // Second wavefunction
    MultiParticleWF waveFunctionTestB(2,2);

    // set it in the state 2*|1>
    waveFunctionTestB(0) = 0;
    waveFunctionTestB(1) = 0;
    waveFunctionTestB(2) = 0;
    waveFunctionTestB(3) = std::complex<double>(2,2);
    BOOST_CHECK(waveFunctionTestB.getDim() == 4);
    BOOST_CHECK(waveFunctionTestB.getN() == 2);
    // Build the density matrix corresponding to the classical mixture  {|0><0|; |1><1|}, with same probability 1/2
    // (0.5, 0, 0, 0)
    // (0, 0, 0, 0)
    // (0, 0, 0, 0)
    // (0, 0, 0, 0.5)
    DensityMatrix densityMatrixB({waveFunctionTestA, waveFunctionTestB}, {0.5, 0.5});
    BOOST_CHECK(densityMatrixB.getDim() == 4);
    BOOST_CHECK(densityMatrixB.getN() == 2);
    CompareComplexNumbers(densityMatrixB.trace(), 1.0);
    for(size_t i = 0; i < 4; i++){
        for(size_t j = 0; j < 4; j++){
            if((i == 0 && j == 0) || (i == 3 && j == 3)){
                CompareComplexNumbers(densityMatrixB(i,j), 0.5);
            } else {
                CompareComplexNumbers(densityMatrixB(i,j), 0.0);
            }
        }
    }

    //Finally test that we get the same result 1/2*|0><0| + 1/2*|1><1| 
    // by using the constructor that takes in input two PureDensityMatrices
    DensityMatrix densityMatrixC({densityMatrixA, PureDensityMatrix(waveFunctionTestB)}, {0.5, 0.5});
    BOOST_CHECK(densityMatrixC.getDim() == 4);
    BOOST_CHECK(densityMatrixC.getN() == 2);
    CompareMatrices(densityMatrixB, densityMatrixC);
}

// Test the function to perform partial trace
BOOST_AUTO_TEST_CASE(density_matrix_partial_trace_test){
    // In the whole test wavefunctions/ density matrices will represent a system of two qubits
    // -----------------------------------------------------
    // ----- TEST 1: perform partial trace of the bell state:
    // |\Psi> =  1/sqrt(2) (|01> - |10>)
    MultiParticleWF waveFunctionTestA(2,2);
    waveFunctionTestA(0) = 0.0;
    waveFunctionTestA(1) = 1.0;
    waveFunctionTestA(2) = -1.0;
    waveFunctionTestA(3) = 0.0;
    waveFunctionTestA.normalize();
    BOOST_CHECK_SMALL(waveFunctionTestA.getNorm() - 1.0, epsilon);

    // The corresponding DensityMatrix is:
    // (0,    0,    0, 0)
    // (0,  0.5, -0.5, 0)
    // (0, -0.5,  0.5, 0)
    // (0,    0,    0, 0)
    DensityMatrix densityMatrixA({waveFunctionTestA}, {1.0});
    // Now perform the partial trace over the subsystem 0 and 1, in this case result is:
    // (0.5,   0)
    // (  0, 0.5)
    // For both subsystems.
    DensityMatrix densityMatrixASubsys0 = densityMatrixA.partialTrace(0);
    DensityMatrix densityMatrixASubsys1 = densityMatrixA.partialTrace(1);
    // Verify that the result is indeed half the identity
    BOOST_CHECK(densityMatrixASubsys0.getN() == 1);
    BOOST_CHECK(densityMatrixASubsys0.getDim() == 2);
    BOOST_CHECK(densityMatrixASubsys1.getN() == 1);
    BOOST_CHECK(densityMatrixASubsys1.getDim() == 2);
    CompareComplexNumbers(densityMatrixASubsys0.trace(), 1.0);
    CompareComplexNumbers(densityMatrixASubsys1.trace(), 1.0);
    for(size_t i = 0; i < 2; i++){
        for(size_t j = 0; j < 2; j++){
            if(i == j){
                CompareComplexNumbers(densityMatrixASubsys0(i,j), 0.5);
                CompareComplexNumbers(densityMatrixASubsys1(i,j), 0.5);
            } else {
                CompareComplexNumbers(densityMatrixASubsys0(i,j), 0.0);
                CompareComplexNumbers(densityMatrixASubsys1(i,j), 0.0);
            }
        }
    }
    // ------------------------------------------------------------
    // ----- TEST 2: perform partial trace of the classical mixture
    // {|01><01| ; |11><11|} with same probability 1/2
    MultiParticleWF waveFunctionTestB_1(2,2);
    waveFunctionTestB_1(0) = 0.0;
    waveFunctionTestB_1(1) = 1.0;
    waveFunctionTestB_1(2) = 0.0;
    waveFunctionTestB_1(3) = 0.0;
    waveFunctionTestB_1.normalize();
    BOOST_CHECK_SMALL(waveFunctionTestB_1.getNorm() - 1.0, epsilon);

    MultiParticleWF waveFunctionTestB_2(2,2);
    waveFunctionTestB_2(0) = 0.0;
    waveFunctionTestB_2(1) = 0.0;
    waveFunctionTestB_2(2) = 0.0;
    waveFunctionTestB_2(3) = 1.0;
    waveFunctionTestB_2.normalize();
    BOOST_CHECK_SMALL(waveFunctionTestB_2.getNorm() - 1.0, epsilon);

    // The corresponding DensityMatrix is:
    // (0,    0, 0,   0)
    // (0,  0.5, 0,   0)
    // (0,    0, 0,   0)
    // (0,    0, 0, 0.5)
    DensityMatrix densityMatrixB({waveFunctionTestB_1, waveFunctionTestB_2}, {0.5, 0.5});
    // Now perform the partial trace over the subsystem 0, the result must be:
    // (0.5,   0)
    // (  0, 0.5)
    DensityMatrix densityMatrixBSubsys0 = densityMatrixB.partialTrace(0);
    BOOST_CHECK(densityMatrixBSubsys0.getN() == 1);
    BOOST_CHECK(densityMatrixBSubsys0.getDim() == 2);
    for(size_t i = 0; i < 2; i++){
        for(size_t j = 0; j < 2; j++){
            if(i == j){
                CompareComplexNumbers(densityMatrixBSubsys0(i,j), 0.5);
            } else {
                CompareComplexNumbers(densityMatrixBSubsys0(i,j), 0.0);
            }
        }
    }
    // Now perform the partial trace over the subsystem 1, the result must be:
    // (0, 0)
    // (0, 1)
    DensityMatrix densityMatrixBSubsys1 = densityMatrixB.partialTrace(1);
    BOOST_CHECK(densityMatrixBSubsys1.getN() == 1);
    BOOST_CHECK(densityMatrixBSubsys1.getDim() == 2);
    for(size_t i = 0; i < 2; i++){
        for(size_t j = 0; j < 2; j++){
            if(i == 1 && j == 1){
                CompareComplexNumbers(densityMatrixBSubsys1(i,j), 1.0);
            } else {
                CompareComplexNumbers(densityMatrixBSubsys1(i,j), 0.0);
            }
        }
    }
}
BOOST_AUTO_TEST_SUITE_END()