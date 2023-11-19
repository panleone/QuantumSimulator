#include "hamiltonian.h"
#include "constants.h"
#include "hermitianMatrix.h"


Hamiltonian::Hamiltonian(std::size_t dim, double leftBound, double rightBound, double mass) 
                            : Hamiltonian(Grid(leftBound, rightBound, dim), mass) {}
Hamiltonian::Hamiltonian(const Grid& grid, double mass) : grid(grid), 
                            H(grid.getDim()), mass(mass)
{
    SetKynetic4Ord();
}

void Hamiltonian::SetKynetic2Ord(){
    const double step = grid.getStep();
    const std::size_t dim = grid.getDim();
    
    const double factor = hslash*hslash/(step*step*2*mass);
    for(size_t i = 0; i < dim-1; i++){
        H(i,i) += 2 * factor;
        H(i+1, i) += -1 * factor;
        H(i, i+1) += -1 * factor;
    }
    H(dim-1, dim-1) += 2 * factor;
}

void Hamiltonian::SetKynetic4Ord(){
    const double step = grid.getStep();
    const std::size_t dim = grid.getDim();

    const double factor = hslash*hslash/(step*step*2*mass); 
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

void Hamiltonian::SetPotential(double (*potFunction)(double, double)){
    const std::size_t dim = grid.getDim();
    for(size_t i = 0; i < dim; i++){
        H(i,i) += potFunction(grid.getNthPoint(i), mass);
    }
}

const HermitianMatrix& Hamiltonian::getHamiltonian() const {
    return H;
}


bool Hamiltonian::solve(char JOBZ){
    return H.findSpectrum(JOBZ);
}

