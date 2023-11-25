#include "hamiltonian.h"
#include "constants.h"
#include "hermitianMatrix.h"
#include "wavefunction.h"
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <ctime>
#include <iostream>
#include <vector>

//TODO: allocate the memory for the hermitian matrix only when the function solve() is called
Hamiltonian::Hamiltonian(const Grid& spaceGrid, double mass, std::function<double(double, double, double)> potFunction) : HermitianMatrix(spaceGrid.getDim()),
                            spaceGrid(spaceGrid), mass(mass), potential(std::move(potFunction))
{}

Hamiltonian::Hamiltonian(const Grid& spaceGrid, const Grid& timeGrid, double mass, std::function<double(double, double, double)> potFunction) :  HermitianMatrix(spaceGrid.getDim()),
                            spaceGrid(spaceGrid), timeGrid(timeGrid), mass(mass), potential(std::move(potFunction))
{}

void Hamiltonian::setKinetic(){
    const double step = spaceGrid.getStep();
    const std::size_t dim = getDim();

    const double factor = hslash*hslash/(step*step*2*mass); 
    auto& H = (*this);
    size_t i = 0;
    for(i = 0; i < dim-2; i++){
        H(i,i) = 5.0/2.0 * factor;
        H(i+1, i) = -4.0/3.0 * factor;
        H(i, i+1) = -4.0/3.0 * factor;
        H(i, i+2) = 1.0/12.0 * factor;
        H(i+2, i) = 1.0/12.0 * factor;
    }
    i = dim-2;
    H(i,i) = 5.0/2.0 * factor;
    H(i+1, i) = -4.0/3.0 * factor;
    H(i, i+1) = -4.0/3.0 * factor;

    i = dim-1;
    H(i,i) = 5.0/2.0 * factor;
}

void Hamiltonian::setPotential(){
    const std::size_t dim = getDim();
    double time = timeGrid.has_value() ? timeGrid->getNthPoint(currentTimeIndex) : 0.0;
    for(size_t i = 0; i < dim; i++){
        (*this)(i,i) += potential(spaceGrid.getNthPoint(i), mass, time);
    }
}

bool Hamiltonian::solve(char JOBZ){
    // TODO: keep track of the last time index at which we solved
    // and reset only if the current hamiltonian time is changed
    reset();
    setKinetic();
    setPotential();
    return findSpectrumAlg2(JOBZ);
}

void Hamiltonian::reset(){
    // If we have have already solved then we need to empty the whole matrix
    if(loadedEigenValues || loadedEigenVectors){
        for(size_t i = 0;i < getDim(); i++){
            for(size_t j = 0; j< getDim(); j++){
                (*this)(i,j) = 0.0;
            }
        }
    }
    loadedEigenValues = false;
    loadedEigenVectors = false;
}

bool Hamiltonian::setTime(size_t newTimeIndex) {
    // If H is time independent bail out early
    if(!timeGrid.has_value()) return false;
    if(newTimeIndex >= timeGrid->getDim()) return false;
    currentTimeIndex = newTimeIndex;
    return true;
}

bool Hamiltonian::incrementTime(){
    return setTime(currentTimeIndex+1);
}

void Hamiltonian::applyTimeEvolutionOperator(WaveFunction& wavefunction) const{
    // Bail out early if no timeGrid has been provided
    if(!timeGrid.has_value()) return;

    if(wavefunction.isFourierTransformed()) wavefunction.FFT();
    const std::complex<double> imUnit = std::complex<double>(0,1);
    const double timeStep = timeGrid->getStep();
    const double t = timeGrid->getNthPoint(this->currentTimeIndex);

    // 1) Apply exp^{-iV*(\Delta t)/2}
    for(size_t i = 0;i < wavefunction.getDim(); i++){
        const double x = wavefunction.getPosition(i);
        wavefunction(i)*= std::exp(-imUnit*potential(x, mass, t)*timeStep/2.0);
    }
    // 2) FFT and apply exp^{-iT*(\Delta t)}
    wavefunction.FFT();
    for(size_t i = 0;i < wavefunction.getDim(); i++){
        const double p = wavefunction.getMomentum(i);
        const double kinEnergy = p*p/(2*mass);
        //std::cerr << kinEnergy << std::endl;
        wavefunction(i)*= std::exp(-imUnit*kinEnergy*timeStep);
    }
    // 3) FFT again and apply exp^{-iV*(\Delta t)/2}
    wavefunction.FFT();
    for(size_t i = 0;i < wavefunction.getDim(); i++){
        const double x = wavefunction.getPosition(i);
        wavefunction(i)*= std::exp(-imUnit*potential(x, mass, t)*timeStep/2.0);
    }
}
// TODO: rewrite this once matrix multiplication is defined 
double Hamiltonian::calculateEnergy(WaveFunction& wavefunction) const{
    std::complex<double> res = 0;
    const double t =  timeGrid.has_value() ? timeGrid->getNthPoint(currentTimeIndex) : 0.0;
    bool isFT = wavefunction.isFourierTransformed();
    if(isFT) wavefunction.FFT();

    // Compute <\psi|V|\psi>
    for(size_t i = 0;i < wavefunction.getDim(); i++){
        const double x = wavefunction.getPosition(i);
        res += std::conj(wavefunction(i))* potential(x, mass, t)*wavefunction(i);
    }
    
    // Compute <\psi|T|\psi> in the fourier space
    wavefunction.FFT();
    for(size_t i = 0;i < wavefunction.getDim(); i++){
        const double p = wavefunction.getMomentum(i);
        res += std::conj(wavefunction(i))* p*p*0.5/(mass)*wavefunction(i);
    }
    // Don't return a different wavefunction from the one inputed
    if (!isFT) wavefunction.FFT();
    return res.real()/pow(wavefunction.getNorm(),2);
}

