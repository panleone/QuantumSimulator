#include <cassert>
#include <complex>
#include <cstddef>
#include <iostream>
#include "hermitianMatrix.h"
#include "linker.h"
#include "matrix.h"

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

    int dim = getDim();

    // Store lower triangle of m1
    char UPLO = 'L';
    
    // will be set to 0 iff we are able to diagonalize the matrix
    int info;

    // work variables needed to the fortran subroutine, they are calculated internally
    int lwork = -1;
    std::vector<std::complex<double>> work(1);
    std::vector<double> rwork(3*dim -1);

    // first call is just to find the right values of work
    zheev_(&JOBZ,&UPLO, &dim, &*matData.begin(),&dim, &*eigenValues.begin(), &*work.begin(), &lwork,&*rwork.begin(),&info);
    // Resize vectors accordingly
    lwork = work.at(0).real();
    work.resize(lwork);
    // Second call is to find eigenvalues/eigenvectors
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

// TODO: and UPLO as inputs
bool HermitianMatrix::findSpectrumAlg2(char JOBZ) {
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

    int dim = getDim();

    // Store lower triangle of m1
    char UPLO = 'L';
    
    // will be set to 0 iff we are able to diagonalize the matrix
    int info;

    // work variables needed to the fortran subroutine, they are calculated internally
    int lwork = -1;
    int lrwork = -1;
    int liwork = -1;
    std::vector<std::complex<double>> work(1);
    std::vector<int> iwork(1);
    std::vector<double> rwork(1);
    // first call is just to find the right values of lwork lrwork and liwork
    zheevd_(&JOBZ,&UPLO, &dim, &*matData.begin(),&dim, &*eigenValues.begin(), &*work.begin(), &lwork,&*rwork.begin(), &lrwork, &*iwork.begin(), &liwork,&info);
    // Resize vectors accordingly
    lwork = work.at(0).real();
    work.resize(lwork);
    lrwork = rwork.at(0);
    rwork.resize(lrwork);
    liwork = iwork.at(0);
    iwork.resize(liwork);
    // Second call is to find eigenvalues/eigenvectors
    zheevd_(&JOBZ,&UPLO, &dim, &*matData.begin(),&dim, &*eigenValues.begin(), &*work.begin(), &lwork,&*rwork.begin(), &lrwork, &*iwork.begin(), &liwork,&info);
    
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

// TODO: remove the first assert in future and load eigenVectors lazily
const ComplexVector HermitianMatrix::getNthEigenVector(size_t N) const {
    assert(loadedEigenVectors && " Hamiltonian hasnt been solved yet");
    assert(N < getDim() && " EigenVector out of bound");
    ComplexVector eigenvector{getDim()};
    for(size_t i = 0; i < getDim(); i++){
        eigenvector(i) = (*this)(i, N);
    }
    return eigenvector;
}

HermitianMatrix ComplexVector::getProjector() const {
    return static_cast<SquareMatrix<std::complex<double>>>(tens_product(*this,getAdjoint()));
}

ComplexMatrix ComplexMatrix::getAdjoint() const {
    ComplexMatrix adjoint{getDimY(), getDimX()};
    for(size_t i = 0; i < getDimX(); i++){
        for(size_t j = 0; j < getDimY(); j++){
        adjoint(j,i) = std::conj((*this)(i,j));
        }
    }
    return adjoint;
}