#include <cassert>
#include <complex>
#include "hermitianMatrix.h"
#include "linker.h"
#include "matrix.h"

std::complex<double> HermitianMatrix::trace() const{
    std::complex<double> trace{0, 0};
    for(std::size_t i = 0; i < getDimX(); i++){
        trace += (*this)(i,i);
    }
    return trace;
}

// TODO: and UPLO as inputs
bool HermitianMatrix::findSpectrum(char JOBZ) {
    // If we alredy fully solved return true
    if(isSolved()) return true;
    // If we don't but we have solved for eigenvalues:
    if(loadedEigenValues){
        if (JOBZ == 'N'){
            return true;
        }else{
            return false;
        }
    }

    int dim = getDimX();

    // Store lower triangle of m1
    char UPLO = 'L';
    
    // will be set to 0 iff we are able to diagonalize the matrix
    int info;

    // work variables needed to the fortran subroutine, we will not use those.
    int lwork = 2*dim -1;
    std::vector<std::complex<double>> work(lwork);
    std::vector<double> rwork(3*dim -1);

    // find eigenvalues
    zheev_(&JOBZ,&UPLO, &dim, &*matData.begin(),&dim, &*eigenValues.begin(), &*work.begin(), &lwork,&*rwork.begin(),&info);
    
    // Update load variables
    bool successful = info == 0;
    if( JOBZ == 'N' ) loadedEigenValues = successful;
    if(JOBZ == 'V') {
        loadedEigenValues = successful;
        loadedEigenVectors = successful;
    }
    return successful;
}

bool HermitianMatrix::isSolved() const {
    return loadedEigenValues && loadedEigenVectors;
}

// TODO: remove the assert in future and load eigenvalues lazily
const std::vector<double>& HermitianMatrix::getEigenValues() const {
    assert(loadedEigenValues && " Hamiltonian hasnt been solved yet");
    return eigenValues;
}
// TODO: remove the assert in future and load eigenVectors lazily
const Matrix<std::complex<double>>& HermitianMatrix::getEigenVectors() const {
    assert(loadedEigenVectors && " Hamiltonian hasnt been solved yet");
    return *this;
}
