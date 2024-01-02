#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <cstddef>
#include <math.h>

#include "test_utils.h"
#include "../grid.h"
#include "../constants.h"
#include "../hamiltonian.h"
#include "../functions.h"


/**
 * @brief Unit tests for the Hamiltonian class
 */
BOOST_AUTO_TEST_SUITE(hamiltonian_tests)
BOOST_AUTO_TEST_CASE(free_hamiltonian_test){
    constexpr double mass = 5.0;
    constexpr double latticeBound = 0.1;
    constexpr size_t dim = 4000;
    Grid grid(0, latticeBound, dim);
    Hamiltonian H(grid, mass, [](double x, double m, double t) { return 0.0;});

    //Test 2: Solve the infinite potential well in the interval (0, 0.1) 
    // and test correctness of eigenvalues/eigenfunctions
    H.solve('V');
    // We are only checking first 1000 eigenvalues since going further would require bigger matrices and bigger run times
    for(size_t i = 0; i < 1000; i++){
        // we are fine if the difference between computed and expected is smaller than 1% of expected result
        BOOST_CHECK_SMALL(H.getEigenValues().at(i) - InfiniteWellEnergy(i+1, latticeBound, mass), InfiniteWellEnergy(i+1, latticeBound, mass)/100);
    }
    // BOOST Check only the first 3 eigenfunctions
    for(size_t i = 0; i < 3; i++){
        for(size_t j = dim/4; j < 3*dim/4; j++){
            double  expected = InfiniteWellWavefunction(i+1, latticeBound, grid.getNthPoint(j));
            int sign = sgn(expected)*sgn(H.getEigenVectors()(j,i).real());
            // Values are in the range [-4.47214, 4.47214]
            // In order to keep the matrix dimension small enough we use a weaker version of epsilon
            // and require that the difference between simulated and expected is smaller than 10^{-2} 
            const double weak_epsilon = epsilon * 100;
            BOOST_CHECK_SMALL(sign*H.getEigenVectors()(j,i).real()/sqrt(grid.getStep()) - expected, weak_epsilon);
        }
    }
}
BOOST_AUTO_TEST_CASE(hamiltonian_test_with_interaction){
    //Same as before but this time we have a harmonic oscillator as interaction
    constexpr double mass = 1.0;
    constexpr double latticeBound = 30;
    constexpr double omega = 1.0;
    constexpr size_t dim = 3000;
    Grid grid(-latticeBound, latticeBound, dim);
    Hamiltonian H(grid, mass, [](double x, double m, double t) { return m*omega*x*x/2;});
    H.solve('V');

    //BOOST check the first 400 eigenvalues
    for(size_t i = 0; i < 400; i ++){
        // we are fine if the difference between computed and expected is smaller than 1% of expected result
        BOOST_CHECK_SMALL(H.getEigenValues().at(i) - HarmonicOscillatorEnergy(i, omega), HarmonicOscillatorEnergy(i, omega)/100);
    }
    // BOOST Check only the first 2 eigenfunctions
    for(size_t i = 0; i < 2; i++){
        for(size_t j = dim/4; j < 3*dim/4; j++){
            double  expected = HarmonicOscillatorWaveFunction(i, grid.getNthPoint(j), mass, omega);
            int sign = sgn(expected)*sgn(H.getEigenVectors()(j,i).real());
            // Values are in the range ~[-1, 1]
            // In order to keep the matrix dimension small enough we use a weaker version of epsilon
            // and require that the difference between simulated and expected is smaller than 0.15
            // Yes it is a big error but to get better results the matrix size would need to be huge 
            const double weak_epsilon = epsilon * 1500;
            BOOST_CHECK_SMALL(sign*H.getEigenVectors()(j,i).real()/sqrt(grid.getStep()) - expected, weak_epsilon);
        }
    }
}
// Create a time dependent harmonic oscillator and test evolution
BOOST_AUTO_TEST_CASE(hamiltonian_evolution_test){
    constexpr double mass = 1.0;
    constexpr double omega = 1.0;

    constexpr double latticeBound = 30;    
    constexpr size_t spaceDim = 3000;

    constexpr double maxTime = 10;
    constexpr double timeDim = 10;
    Grid spaceGrid(-latticeBound, latticeBound, spaceDim);
    Grid timeGrid(0, maxTime, timeDim);

    // This is a free hamiltonian at t=0 and becomes a normal harmonic oscillator at t=tMax
    Hamiltonian H(spaceGrid, timeGrid, mass, [](double x, double m, double t) { return t/maxTime*m*omega*x*x/2;});

    // Solve and verify that we found the expected eigenvalues for free hamiltonian
    BOOST_CHECK(H.solve('N'));
    for(size_t i = 0; i < 400; i ++){
        // we are fine if the difference between computed and expected is smaller than 1% of expected result
        // Energies are in the range [0,220]
        BOOST_CHECK_SMALL(H.getEigenValues().at(i) - InfiniteWellEnergy(i+1, 2*latticeBound, mass), InfiniteWellEnergy(i+1, 2*latticeBound, mass)/100);
    }
    // Now evolve it, we can do it for at most timeDim-1 times:
    for(size_t i = 0; i < timeDim-1; i ++){
        BOOST_CHECK(H.incrementTime());
    }
    BOOST_CHECK(!H.incrementTime());
    // Solve again and verify that we get harmonic hoscillator solutions
    BOOST_CHECK(H.solve('N'));
    for(size_t i = 0; i < 400; i ++){
        BOOST_CHECK_SMALL(H.getEigenValues().at(i) - HarmonicOscillatorEnergy(i, omega), HarmonicOscillatorEnergy(i, omega)/100);
    }
}
BOOST_AUTO_TEST_SUITE_END()