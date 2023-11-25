// This file contains useful functions that are used in different parts of the project
#ifndef FUNCTIONS 
#define FUNCTIONS
#include <cstddef>
#include <math.h>
#include "constants.h"

/**
 * @brief computes n-th energy level for an infinite well  of length L
 * for a particle of mass m
 */
inline double InfiniteWellEnergy(size_t n, double L, double m){
    return n*n*hslash*hslash*M_PI*M_PI/(2*m*L*L);
}

/**
 * @brief Computes the n-th eigenfunction for an infinite well in the range [0,L]
 * evaluated at point x
 */
inline double InfiniteWellWavefunction(size_t n, double L, double x){
    return sqrt(2/L)*sin(x*n*M_PI/L);
}

/**
 * @brief computes n-th energy level for a harmonic oscillator of pulse omega
 * for a particle of mass m 
 */
inline double HarmonicOscillatorEnergy(size_t n, double omega){
    return hslash*omega*(n +1/2.0);
}

/**
 * @brief Computes the n-th theoretical eigenvector for the harmonic oscillator at point x
 * 
 * @param n - n-th eigenstate \ket{\psi_{n}}
 * @param x - coordinate at which we want to evaluate \psi_{n}
 * @return double \psi_{n}(x)
 */
inline double HarmonicOscillatorWaveFunction(int n, double x, double mass, double omega){
    const double param = mass*omega/hslash;
    const double factor1 = pow(sqrt(pow(2, n)*tgamma(n+1)), -0.5);
    const double factor2 = exp(-param*x*x/2);
    const double factor3 = pow(param/M_PI, 0.25);
    return factor1 * factor2 * factor3 * std::hermite(n, x*sqrt(param));
}

#endif // FUNCTIONS