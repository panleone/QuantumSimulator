#include <cmath>
#include <cstddef>
#include <iostream>
#include <math.h>
#include "constants.h"
#include "hamiltonian.h"
#include "timer.h"
#include <fstream>

/**
 * @brief Computes the n-th theoretical eigenvector for the harmonic oscillator at point x
 * 
 * @param n - n-th eigenstate \ket{\psi_{n}}
 * @param x - coordinate at which we want to evaluate \psi_{n}
 * @return double \psi_{n}(x)
 */
double teoSol(int n, double x, double mass, double omega){
    const double param = mass*omega/hslash;
    const double factor1 = pow(sqrt(pow(2, n)*tgamma(n+1)), -0.5);
    const double factor2 = exp(-param*x*x/2);
    const double factor3 = pow(param/M_PI, 0.25);
    return factor1 * factor2 * factor3 * std::hermite(n, x*sqrt(param));
}

template<typename T>
/**
 * @brief Parse the input argv inside parsedVar
 * @return true iff the operation was succesful
 */
bool parseInput(T& parsedVar, const std::string& argv){
    std::istringstream ss(argv);
    if (!(ss >> parsedVar)) {
        std::cerr << "Invalid arg1: " << argv << std::endl;
        return false;
    } else if (!ss.eof()) {
        std::cerr << "Trailing characters after number: " << argv << std::endl;
        return false;
    }
    return true;
}

/**
 * @brief Solves the Quantum Harmonic Oscillator and outputs eigenvectors and eigenvalues.
 */
int main(int argc, char **argv){
    
    // Correctly parse the input
    if ( argc != 3){
        std::cerr << "Run the program with only two argument!\n\
-first argument  arg1 (positive double) the interval extremes [-arg1, arg1] on which the shrodinger equation is solved\n\
-second argument arg2 (uint) the discretization dimension: the interval will be discretized in arg2 steps" << std::endl;
        return 0;
    }

    double latticeBound;
    std::size_t dim;
    if(!parseInput(latticeBound, argv[1]) || !parseInput(dim, argv[2])) return 0;
    if(dim<2){
        std::cerr << "arg2 must be at LEAST 2" <<std::endl;
        return 0;
    }
    if(latticeBound < 0){
        std::cerr << "arg1 MUST BE positive" << std::endl;
    }

    // Parameters of our problem
    constexpr double mass = 1.0;
    constexpr double omega = 1.0;
    Timer timer;
    // Grid on which we are solving the shrodinger equation
    Grid grid(-latticeBound, latticeBound, dim);

    // Create the Hamiltonian and set potential with a lambda
    Hamiltonian H1(grid, mass);
    H1.setPotential( [](double x, double m) { return omega*m*x*x/2;});

    std::cerr << "Time elapsed to set everything up: " << timer.elapsed() << " seconds\n";
    H1.solve('V');
    std::cerr << "Time elapsed to solve: " << timer.elapsed() << " seconds\n";

    // Save data on file
    for(size_t i = 0; i < dim; i ++){
        std::cout << H1.getHamiltonian().getEigenValues().at(i) << " " ;
    }
    std::cout << "\n";
    // Save eigenfunctions on file
    int selected;    
    while(true){
        std::cerr << "Which eigenfunction you want to output? (insert -1 to terminate)" << std::endl;
        std::cin >> selected;
        if(selected < 0 || selected >= dim){
            break;
        }
        for(size_t i = 0; i < dim; i ++){
            std::cout << H1.getHamiltonian().getEigenVectors()(i, selected) << " " ;
        }
        std::cout << "\n";
    }
return 0;
}