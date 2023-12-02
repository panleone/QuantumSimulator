#include <cassert>
#include <complex>
#include <cstddef>
#include <iostream>
#include <vector>

#include "functions.h"
#include "hamiltonian.h"
#include "wavefunction.h"

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
    // Discretize the time by dividing [0, maxT] in timeDim points
    constexpr double maxT = 10*3;
    constexpr size_t timeDim = 5000*3;
    Grid timeGrid = Grid(0,maxT,timeDim);

    // Discretize the space by dividing [-spaceGridBound, spaceGridBound] in spaceDim points
    constexpr size_t spaceDim = 20000;
    constexpr double maxSpace = 40;
    Grid spaceGrid = Grid(-maxSpace,maxSpace,spaceDim);

    //We are considering a moving potential in the range [0, potFinalPosition]
    constexpr double potFinalPosition = 10;
    // particle and potential parameters
    constexpr double mass = 1.0;
    constexpr double omega = 1.0;

    Hamiltonian H(Grid(-maxSpace, maxSpace, 2), timeGrid, mass, 
                    [](double x, double m, double t) { return m*omega*pow(x -potFinalPosition*t/maxT,2)/2;});
    
    WaveFunction wavefunction(spaceGrid);
    // Vector where the average positions of the W.F. are saved
    std::vector<double> posAverage(timeDim-1);
    std::vector<double> posVariance(timeDim-1);
    // Number of snapshoots of the wave function we are taking
    constexpr size_t nSnapshoots = 600;
    Matrix<double> wavefunctionSnapshoots(nSnapshoots, spaceDim);

    // At time t = 0 our wavefunction is in the 0-th eigenvalue of H(t=0)
    for(size_t i = 0; i < spaceDim; i++){
        wavefunction(i) = HarmonicOscillatorWaveFunction(0, spaceGrid.getNthPoint(i), mass, omega);
    }
    wavefunction/= wavefunction.getNorm();

    // Evolve the system
    for(size_t i = 0; i < timeDim -1;i++){
        H.applyTimeEvolutionOperator(wavefunction);
        assert(H.incrementTime());
        posAverage.at(i) = wavefunction.positionAverage();
        posVariance.at(i) = wavefunction.positionVariance();

        // Take a snapshoot of the wave function every (timeDim/nSnapshoots) time intervals
        if(i % (timeDim/nSnapshoots) == 0){
            size_t currSnapshot = i/(timeDim/nSnapshoots);
            std::cerr << "Snapshoot taken:" << " " <<currSnapshot << " of " << nSnapshoots<< std::endl;
            for(size_t j = 0; j < spaceDim; j++){
                wavefunctionSnapshoots(currSnapshot, j) = std::norm(wavefunction(j));
            }
        }
    }
    // cout the average positions and variance
    for(auto& x : posAverage){
        std::cout << x << " ";
    }
    std::cout << "\n";
    for(auto& x : posVariance){
        std::cout << x << " ";
    }
    std::cout << "\n";
    std::cout << wavefunctionSnapshoots << std::endl;
return 0;
}