#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <cstddef>

#include "test_utils.h"
#include "../grid.h"
#include "../constants.h"
#include "../hamiltonian.h"
#include "math.h"

/**
 * @brief computes n-th energy level for an infinite well in the range [0, L]
 * for a particle of mass m
 */
double InfiniteWellEnergy(size_t n, double L, double m){
    return n*n*hslash*hslash*M_PI*M_PI/(2*m*L*L);
}

/**
 * @brief Computes the n-th eigenfunction for an infinite well in the range [0,L]
 * evaluated at point x
 */
double InfiniteWellWavefunction(size_t n, double L, double x){
    return sqrt(2/L)*sin(x*n*M_PI/L);
}

/**
 * @brief Unit tests for the Hamiltonian class
 */
BOOST_AUTO_TEST_SUITE(hamiltonian_tests)
BOOST_AUTO_TEST_CASE(hamiltonian_test){
    constexpr double mass = 5.0;
    constexpr double latticeBound = 0.1;
    constexpr size_t dim = 8000;
    Grid grid(0, latticeBound, dim);
    Hamiltonian H(grid, mass);

    //Test 1: Test the correctness of the kinetic term
    std::complex<double> expected;
    const double factor = hslash*hslash/(2*mass*grid.getStep()*grid.getStep());
    for(size_t i = 0; i < dim; i++){
        for(size_t j = 0; j < dim; j++){
            if(i == j){
                expected = 5.0/2.0 * factor;
            }else if(i == j+1 || j == i+1){
                expected = -4.0/3.0 * factor;
            }else if(i == j+2 || j == i+2){
                expected = 1.0/12.0 * factor;
            }else{
                expected = 0.0;
            }
            CompareComplexNumbers(H.getHamiltonian()(i,j), expected);
        }
    }
    //Test 2: Solve the infinite potential well in the interval (0, 0.1) 
    // and test correctness of eigenvalues/eigenfunctions
    H.solve('V');
    // We are only checking first 2000 eigenvalues since going further would require bigger matrices and bigger run times
    for(size_t i = 0; i < 2000; i++){
        // we are fine if the difference between computed and expected is smaller than 1% of expected result
        BOOST_CHECK_SMALL(H.getHamiltonian().getEigenValues().at(i) - InfiniteWellEnergy(i+1, latticeBound, mass), InfiniteWellEnergy(i+1, latticeBound, mass)/100);
    }
    // BOOST Check only the first 10 eigenfunctions
    for(size_t i = 0; i < 10; i++){
        for(size_t j = 1000; j < 4000; j++){
            double  expected = InfiniteWellWavefunction(i+1, latticeBound, grid.getNthPoint(j));
            int sign = sgn(expected)*sgn(H.getHamiltonian().getEigenVectors()(j,i).real());
            // Values are in the range [-4.47214, 4.47214]
            // In order to keep the matrix dimension small enough we use a weaker version of epsilon
            // and require that the difference between simulated and expected is smaller than 10^{-2} 
            const double weak_epsilon = epsilon * 100;
            BOOST_CHECK_SMALL(sign*H.getHamiltonian().getEigenVectors()(j,i).real()/sqrt(grid.getStep()) - expected, weak_epsilon);
        }
    }
}
BOOST_AUTO_TEST_SUITE_END()