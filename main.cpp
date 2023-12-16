#include <cassert>
#include <complex>
#include <cstddef>
#include <iostream>
#include <vector>
#include "timer.h"
#include "hermitianMatrix.h"

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

int main(int argc, char **argv){
    // Correctly parse the input
    if ( argc != 3){
        std::cerr << "Run the program with only two argument!\n\
-first argument  arg1 (positive integer) the number of particles N that will be simulated, must be N>=2\n\
-second argument arg2 (double) strength \\lambda of the external magnetic field " << std::endl;
        return 0;
    }

    // Number of particles
    size_t N = 0;
    // Strength of external field
    double lambda = 0.0;
    if(!parseInput(N, argv[1]) || !parseInput(lambda, argv[2])) return 0;
    if(N<2){
        std::cerr << "arg1 must be at LEAST 2" <<std::endl;
        return 0;
    }

    // Build identity matrix
    HermitianMatrix identity(2);
    identity(0,0) = identity(1,1) = 1.0;

    // Build sigma x
    HermitianMatrix sigmaX(2);
    sigmaX(0,1) = sigmaX(1,0) = 1.0;

    // Build sigma z
    HermitianMatrix sigmaZ(2);
    sigmaZ(0,0) = 1.0;
    sigmaZ(1,1) = -1.0;

    // Build the Hamiltonian, it has size 2^N
    HermitianMatrix H(pow(2, N));
    Timer timer{};

    // Add the Interaction with external magnetic field:
    for(size_t i = 0; i < N; i++){
        HermitianMatrix interaction = i == 0 ? sigmaZ : identity;
        for(size_t j = 1; j < N; j++){
            interaction = tens_product(interaction, ((j == i) ? sigmaZ : identity));
        }
        H += interaction* static_cast<std::complex<double>>(lambda);
    }
    // Add the Interaction among nearest neighbor spins:
    for(size_t i = 0; i < N-1; i++){
        HermitianMatrix interaction = i == 0 ? sigmaX : identity;
        for(size_t j = 1; j < N; j++){    
            interaction = tens_product(interaction, ((j == i || j == i+1) ? sigmaX : identity)); 
        }
        H += interaction;
    }
    std::cerr << "Hamiltonian has been built, diagonalizing:" << std::endl;
    // Find and print the spectrum
    H.findSpectrumAlg2('N');
    for(auto& E : H.getEigenValues()){
        std::cout << E << " , ";
    }
    std::cout << "\n";
    std::cerr << "Total time elapsed: " << timer.elapsed() << std::endl;
    return 0;
}