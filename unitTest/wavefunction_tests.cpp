#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <complex>
#include <math.h>

#include "test_utils.h"
#include "../grid.h"
#include "../constants.h"
#include "../wavefunction.h"
#include "../hamiltonian.h"
#include "../functions.h"

/**
 * @brief Unit tests for the Wavefunction class
 */
BOOST_AUTO_TEST_SUITE(wavefunction_tests)

// Test that we are storing the right momentum in the wafe function momentum grid
// and that we are doing fourie transform properly
BOOST_AUTO_TEST_CASE(fourierTransform_test){
    // Create a grid
    constexpr double leftBound = -40.52;
    constexpr double rightBound = 3100.45;
    std::size_t dim = 10;
    Grid grid(leftBound, rightBound, dim);
    // On which we define a (normalized) wave function
    WaveFunction testWF(grid);
    for(int i = 0; i< dim; i++){
        testWF(i) = static_cast<std::complex<double>>(i);
    }
    testWF /= testWF.getNorm();
    // Save in memory a copy
    const WaveFunction testWFCopy = testWF;
    // Check that initial setup is good
    BOOST_CHECK_SMALL(testWF.getPosition(0) - grid.getLeftBound(), epsilon);
    BOOST_CHECK_SMALL(testWF.getPosition(dim - 1) - grid.getRightBound(), epsilon);
    BOOST_CHECK(!testWF.isFourierTransformed());
    BOOST_CHECK_SMALL(testWF.getNorm() - 1.0 , epsilon);
    BOOST_CHECK(testWF.getDim() == dim);

    // Apply FFT
    testWF.FFT();
    BOOST_CHECK_SMALL(testWF.getNorm() - 1.0 , epsilon);
    BOOST_CHECK(testWF.isFourierTransformed());

    // Verify that we obtain the same result by manually doing the FFT
    for(int i = 0; i<dim; i++){
        std::complex<double> res = 0;
        // Momentum of the i-th component of the WF
        double pi = testWF.getMomentum(i);
        for(int j = 0; j<dim; j++){
            // Coordinate (shifted) of the j-th component of the WF
            double xj = testWF.getPosition(j) - testWF.getPosition(0);
            res += testWFCopy(j)*std::exp(-std::complex<double>(0,1)*pi*xj)/sqrt(dim);
        }
        CompareComplexNumbers(res, testWF(i));
    }

    // Test the momenta given by the function getMomentum() are correct by doing manually inverse ft
    // snd verifying that we get the original wave function in coordinate space
    for(int i = 0; i<dim; i++){
        std::complex<double> res = 0;
        double xi = testWF.getPosition(i) - testWF.getPosition(0);
        for(int j = 0; j<dim; j++){
            double pj = testWF.getMomentum(j);
            res += testWF(j)*std::exp(std::complex<double>(0,1)*pj*xi)/sqrt(dim);
        }
        // we should have the original wave function
        CompareComplexNumbers(res, testWFCopy(i));
    }

    // Finally make the lib to a inverse FFT and verify that we get back the original vector
    testWF.FFT();
    BOOST_CHECK_SMALL(testWF.getNorm() - 1.0 , epsilon);
    BOOST_CHECK(!testWF.isFourierTransformed());
    CompareMatrices(testWF, testWFCopy);
}
// Test average and variance functions
BOOST_AUTO_TEST_CASE(average_variance_test){

    constexpr double latticeBound = 6;
    constexpr double latticeDim =  4;
    // Grid [0.0, 2.0, 4.0, 6.0]
    Grid spaceGrid(0, latticeBound, latticeDim);
    WaveFunction testWF(spaceGrid);
    testWF(0) = 0;
    testWF(1) = 0;
    testWF(2) = 100;
    testWF(3) = 0;

    BOOST_CHECK_SMALL(testWF.positionAverage() - spaceGrid.getNthPoint(2), epsilon);
    BOOST_CHECK_SMALL(testWF.positionVariance(), epsilon);
    BOOST_CHECK_SMALL(testWF.getNorm() - 100, epsilon);

    testWF(0) = 100;
    BOOST_CHECK_SMALL(testWF.positionAverage() - spaceGrid.getNthPoint(1), epsilon);
    BOOST_CHECK_SMALL(testWF.positionVariance() - 4, epsilon);
    BOOST_CHECK_SMALL(testWF.getNorm() - sqrt(2*pow(100,2)), epsilon);

    testWF(3) = 50;
    BOOST_CHECK_SMALL(testWF.positionAverage() - 2.4444, epsilon);
    BOOST_CHECK_SMALL(testWF.positionVariance() - 5.1358, epsilon);
    BOOST_CHECK_SMALL(testWF.getNorm() - sqrt(2*pow(100,2)+ pow(50,2)), epsilon);
}

// test that the norm is an invariance under time evolution
BOOST_AUTO_TEST_CASE(invariance_of_norm){
    // Discretize the time by dividing [0, maxT] in timeDim points
    constexpr double maxT = 1;
    constexpr size_t timeDim = 5000;
    Grid timeGrid = Grid(0,maxT,timeDim);

    // Discretize the space by dividing [-spaceGridBound, spaceGridBound] in spaceDim points
    constexpr size_t spaceDim = 5000;
    constexpr double spaceGridBound = 10;
    Grid spaceGrid = Grid(-spaceGridBound,spaceGridBound,spaceDim);

    constexpr double mass = 1.0;

    Hamiltonian H(spaceGrid, timeGrid, mass, [](double x, double m, double t) { return t*x;});
    WaveFunction wavefunction(spaceGrid);
    // Pick a generic initial position, like the eigenstate of harmonic oscillator
    for(int i = 0; i < spaceDim; i++){
        wavefunction(i) = HarmonicOscillatorWaveFunction(2, spaceGrid.getNthPoint(i), mass, 1.0);
    }
    double initialNorm = wavefunction.getNorm();

    // Evolve the system
    for(int i = 0;i < timeDim -1;i++){
        H.applyTimeEvolutionOperator(wavefunction);
        BOOST_CHECK(H.incrementTime());
        // initialNorm is 26.58
        BOOST_CHECK_SMALL(wavefunction.getNorm() - initialNorm, epsilon);
    }
}

// test that we get the correct average energy <\psi|H|\psi>
// The calculation is done internally in Fourier space so this is mostly a test for fourier transform
BOOST_AUTO_TEST_CASE(average_energy_test){
    // Prepare a harmonic oscillator
    constexpr size_t spaceDim = 5000;
    constexpr double spaceGridBound = 10;
    Grid spaceGrid = Grid(-spaceGridBound,spaceGridBound,spaceDim);
    constexpr double mass = 1.0;
    constexpr double omega = 1.0;
    Hamiltonian H(Grid(-spaceGridBound, spaceGridBound, 2), mass, [](double x, double m, double t) { return omega*m*pow(x,2)/2;});
    WaveFunction wavefunction(spaceGrid);
    // Pick a wavefunction in the first eigenstate
    for(int i = 0; i < spaceDim; i++){
        wavefunction(i) = HarmonicOscillatorWaveFunction(0, spaceGrid.getNthPoint(i), mass, omega);
    }
    // Boost check that the average energy is 0.5
    BOOST_CHECK_SMALL(H.calculateEnergy(wavefunction) - HarmonicOscillatorEnergy(0, omega), epsilon);
    // Re-do with n=5
    for(int i = 0; i < spaceDim; i++){
        wavefunction(i) = HarmonicOscillatorWaveFunction(5, spaceGrid.getNthPoint(i), mass, omega);
    }
    // Average should be 5.5
    BOOST_CHECK_SMALL(H.calculateEnergy(wavefunction) - HarmonicOscillatorEnergy(5, omega), epsilon);
}
BOOST_AUTO_TEST_SUITE_END()