#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <vector>
#include <math.h>
#include "densityMatrix.h"
#include "matrix.h"
#include "multiParticleWF.h"
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
// Returns the MxM identity matrix
HermitianMatrix Identity(size_t M){
    HermitianMatrix identity{M};
    for(size_t i = 0; i<M; i++){
        identity(i,i) = 1.0;
    }
    return identity;
}

static HermitianMatrix sigmaX(2);
static HermitianMatrix sigmaZ(2);

HermitianMatrix buildHamiltonian(size_t N, double lambda){
    auto identity = Identity(2);
    // Build the Hamiltonian, it has size 2^N
    HermitianMatrix H(pow(2, N));

    // Add the Interaction with external magnetic field:
    HermitianMatrix interaction = identity;
    for(size_t i = 0; i < N; i++){
        interaction = i == 0 ? sigmaZ : identity;
        for(size_t j = 1; j < N; j++){
            interaction = tens_product(interaction, ((j == i) ? sigmaZ : identity));
        }
        interaction*= static_cast<std::complex<double>>(lambda);
        H += interaction;
    }
    // Add the Interaction among nearest neighbor spins:
    for(size_t i = 0; i < N-1; i++){
        interaction = i == 0 ? sigmaX : identity;
        for(size_t j = 1; j < N; j++){    
            interaction = tens_product(interaction, ((j == i || j == i+1) ? sigmaX : identity)); 
        }
        H += interaction;
    }
    return H;
}

/**
 * @brief Computes the ground state of the quantum ising model 
 using the real space RG algorithm given the number of particles N and the field strength \lambda
 * @param N The number of particles
 * @param lambda The strength of the field
 */
void realSpaceRG(size_t N, double lambda){
    // Max iterations that we do before giving up
    constexpr size_t maxIterations = 100;
    // last value for energy density that we found with RG
    double oldEnergyDens = -9999.0;
    // current energy value density found with the RG algorithm
    double currEnergyDens = -9999.0;
    // energy threshold: we have convergence when |currEnergy - oldEnergy|/N < threshold
    double threshold = 0.0001;

    // -----STEP 0: BUILD INITIAL QUANTITIES FOR A,B,H-----
    ComplexMatrix A{tens_product(Identity(pow(2,N-1)),sigmaX)};
    ComplexMatrix B{tens_product(sigmaX, Identity(pow(2,N-1)))};
    HermitianMatrix H = buildHamiltonian(N, lambda);
    HermitianMatrix Hcopy{H};

    for(size_t k = 0; k < maxIterations; k++){
        // -----STEP 1: DOUBLE THE SIZE-----
        H = static_cast<SquareMatrix<std::complex<double>>>(tens_product(H, Identity(pow(2,N))) 
                                                        + tens_product(Identity(pow(2,N)), H) + tens_product(A,B));
        Hcopy = H;
        // -----STEP 2: DIAGONALIZE AND CHECK IF THE ALGORITHM CONVERGED-----  
        H.findSpectrumAlg2('V');
        oldEnergyDens = currEnergyDens;
        currEnergyDens = H.getEigenValues()[0]/(2*N);
        if(abs(oldEnergyDens-currEnergyDens) < threshold){
            std::cerr << "Algorithm converged in: " << k+1 << " iterations!" << std::endl;
            std::cout << "VICTORY! " << currEnergyDens << std::endl; 
            return;
        }

        // -----STEP 3: BUILD THE PROJECTOR ON THE SMALLEST N EIGENVALUES----- 
        ComplexMatrix P{static_cast<size_t>(pow(2,2*N)),static_cast<size_t>(pow(2,N))};
        for(size_t i = 0; i < pow(2,2*N); i++){
            for(size_t j=0; j < pow(2,N); j++){
                P(i,j) += H.getEigenVectors()(i,j);
            }
        }
        // -----STEP 4: UPDATE H,A,B-----
        H = static_cast<SquareMatrix<std::complex<double>>>(matrix_product(P.getAdjoint(), matrix_product(Hcopy, P)));
        H /= 2.0;
        A = matrix_product(P.getAdjoint(),
                matrix_product(
                    tens_product(Identity(pow(2,N)),A)
                    ,P)
            );
        A /= sqrt(2.0);

        B = matrix_product(P.getAdjoint(),
                matrix_product(
                    tens_product(B, Identity(pow(2,N)))
                    ,P)
            );
        B /= sqrt(2.0);
    }
    std::cerr << "Algorithm failed! Energies did not converge." << std::endl;
}

/**
 * @brief Computes the ground state of the quantum ising model 
 using the infinite density matrix RG algorithm given the field strength \lambda
 * @param lambda The strength of the field
 */
void infiniteDensityMatrixRG(std::complex<double> lambda){
    // Max iterations that we do before giving up
    constexpr size_t maxIterations = 500;
    // last value for energy density that we found with RG
    double oldEnergyDens = -9999.0;
    // current energy value density found with the RG algorithm
    double currEnergyDens = -9999.0;
    // energy threshold: we have convergence when |currEnergy - oldEnergy|/N < threshold
    double threshold = 0.0001;
    //Let's start with a system containing two particles.
    // The left one is the left system, the right one is the right system.

    // -----STEP 0: BUILD THE INITIAL HAMILTONIAN WITH ONLY 2 PARTICLES-----
    // Hamiltonian for the Left subsystem
    HermitianMatrix HL{sigmaZ*lambda};
    // Hamiltonian for the Right subsystem
    HermitianMatrix HR{sigmaZ*lambda};

    HermitianMatrix BL {sigmaX};
    HermitianMatrix BR {sigmaX};
    // Interaction among Left and Right subsystems
    HermitianMatrix HLR {tens_product(tens_product(Identity(2),sigmaX),tens_product(sigmaX,Identity(2)))};
    for(size_t k = 0; k < maxIterations; k++){

        // -----STEP 1: EXPAND LEFT AND RIGHT SUBSYSTEMS BY ADDING ONE PARTICLE-----
        HL = static_cast<SquareMatrix<std::complex<double>>>(tens_product(HL,Identity(2)) + tens_product(Identity(2), sigmaZ*lambda) + tens_product(BL,sigmaX));
        HR = static_cast<SquareMatrix<std::complex<double>>>(tens_product(Identity(2), HR) + tens_product(sigmaZ*lambda, Identity(2)) + tens_product(sigmaX,BR));

        BL = tens_product(Identity(2),sigmaX);
        BR = tens_product(sigmaX, Identity(2));
        
        // -----STEP 2: BUILD THE TOTAL HAMILTONIAN AND DIAGONALIZE IT-----
        HermitianMatrix Htot{tens_product(HL, Identity(pow(2,2))) + tens_product(Identity(pow(2,2)), HR) + HLR};
        Htot.findSpectrumAlg2('V');
        oldEnergyDens = currEnergyDens;
        size_t totParticles = (k+1+1)*2;
        currEnergyDens = Htot.getEigenValues().at(0)/(totParticles);
        if(abs(oldEnergyDens-currEnergyDens) < threshold){
            std::cerr << "Algorithm converged in: " << k+1 << " iterations!" << std::endl;
            std::cout << "VICTORY! " << currEnergyDens << std::endl; 
            return;
        }
        // -----STEP 3: BUILD THE DENSITY MATRIX AND TRACE OUT THE RIGHT SUBSYSTEM-----
        MultiParticleWF waveFunction{4,2};
        auto groundState = Htot.getNthEigenVector(0);
        for(size_t i = 0; i < pow(2,4); i++){
            waveFunction(i) = groundState(i);
        }
        DensityMatrix gsDensityMatrix{PureDensityMatrix{waveFunction}};
        // -----STEP 4: DIAGONALIZE the reduced density matrix and build the Projector-----
        // Trace out the right subsystem
        DensityMatrix reducedMatrix = gsDensityMatrix.partialTrace(0);
        reducedMatrix.findSpectrumAlg2('V');
        ComplexMatrix P{4,2};
            for(size_t i = 0; i < 4; i++){
                for(size_t j=0; j < 2; j++){
                    P(i,j) += reducedMatrix.getEigenVectors()(i,2+j);
                }
            }
        // NB: there is surely a way to map the left subsystem to the right one by symmetry arguments
        // I could not figure out how it is done,
        // So I am just computing it in the hard way:
        // i.e. by tracing out the left subsystem and calculating the right projector 
        DensityMatrix reducedMatrixR = gsDensityMatrix.partialTrace(1);
        reducedMatrixR.findSpectrumAlg2('V');
        ComplexMatrix PR{4,2};
            for(size_t i = 0; i < 4; i++){
                for(size_t j=0; j < 2; j++){
                    PR(i,j) += reducedMatrixR.getEigenVectors()(i,2+j);
                }
            }
        
        // -----STEP 5: UPDATE EVERYTHING-----
        HL = static_cast<SquareMatrix<std::complex<double>>>(matrix_product(P.getAdjoint(), matrix_product(HL,P)));
        BL = static_cast<SquareMatrix<std::complex<double>>>(matrix_product(P.getAdjoint(), matrix_product(BL, P)));
        
        // See the comment above in step 4 about computing the new right subsystem
        HR = static_cast<SquareMatrix<std::complex<double>>>(matrix_product(PR.getAdjoint(), matrix_product(HR,PR)));
        BR = static_cast<SquareMatrix<std::complex<double>>>(matrix_product(PR.getAdjoint(), matrix_product(BR, PR)));
    }
    std::cerr << "Algorithm failed! Energies did not converge." << std::endl;
}

int main(int argc, char **argv){
    // Correctly parse the input
    std::cerr << "This programs computes the ground state of the Quantum Ising Model using either the real space RG or the infinite density matrix RG" << std::endl;  
    if ( argc != 4){
        std::cerr << "Run the program with only two argument!\n\
-first argument  arg1 (positive integer) the number of particles N that will be simulated, must be N>=2\n\
-second argument arg2 (double) strength \\lambda of the external magnetic field\n\
-third argument arg3 algorithm mode: insert 0 to run the real space RG, else it will run the infinite density matrix RG" << std::endl;
        return 0;
    }

    // Number of particles
    size_t N = 0;
    // Strength of external field
    double lambda = 0.0;
    // Algorithm mode: 0 to use real space RG, other values will run the infinite density matrix RG.
    int mode = 0;
    if(!parseInput(N, argv[1]) || !parseInput(lambda, argv[2]) || !parseInput(mode, argv[3])) return 0;
    if(N<2){
        std::cerr << "arg1 must be at LEAST 2" <<std::endl;
        return 0;
    }
    // Build sigma x
    sigmaX(0,1) = sigmaX(1,0) = 1.0;
    // Build sigma z
    sigmaZ(0,0) = 1.0;
    sigmaZ(1,1) = -1.0;
    if(mode == 0){
        realSpaceRG(N, lambda);
    }else{
        infiniteDensityMatrixRG(lambda);
    }
    return 0;
}   