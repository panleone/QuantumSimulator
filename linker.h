#ifndef LINKER
#define LINKER
#include <complex>

// This file contains the linker functions between LAPACK (fortran) and C++
// See https://www.netlib.org/lapack/#_documentation for more info.

// Computes LU factorization of A
extern "C" void zgetrf_(int* dim1, int* dim2, std::complex<double>* A, int* lda, int* ipiv, int* info);
// Solves the linear system A*x = B
extern "C" void zgetrs_(char *TRANS, int *N, int *NRHS, std::complex<double> *A, int *LDA, int *IPIV, std::complex<double> *B, int *LDB, int *INFO );
// Computes the eigenvalues of the hermitian matrix A
extern "C" void zheev_( char* jobz, char* uplo, int* n, std::complex<double>* A, int* lda,
                double* w, std::complex<double>* work, int* lwork, double* rwork, int* info );
#endif //LINKER
