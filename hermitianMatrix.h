#ifndef HERMITIAN
#define HERMITIAN

#include "matrix.h"
#include <complex>
#include <cstddef>
#include <vector>

/**
 * TODO: optimize, since the matrix is hermitian we can get rid of almost half of the memory
 * @brief class for Hermitian Matrices with useful functions
 */
class HermitianMatrix : public SquareMatrix<std::complex<double>>
{
    protected:
        //True if we eigenvalues of the hamiltonian have been computed,
        bool loadedEigenValues = false;
        // Eigen values of the matrix, they are loaded lazily or by calling solve()
        std::vector<double> eigenValues;

        //True if we eigenvalues of the hamiltonian have been computed.
        bool loadedEigenVectors = false;
    public:
        explicit HermitianMatrix(std::size_t dim) : SquareMatrix<std::complex<double>>(dim), eigenValues(dim)
        {};
        /**
         * @brief Find eigenvalues and eigenvectors of the hermitian matrix:
         * AFTER CALLING THIS FUNCTION THE MATRIX CONTENT IS DESTROYED
         * in particular if eigenvalues are computed it's content is replaced
         * with the eigenvector matrix
         * @param egv - vector in which eigenvalues will be stored
         * @param JOBZ - char: if set to 'N' will compute eigenvalue only,
         * if set to 'V' will compute also eigenvectors  
         * @return true iff the opeartion was successful
         */
        bool findSpectrum(char JOBZ);
        // Same as findSpectrum but uses a different algorithm
        bool findSpectrumAlg2(char JOBZ);
        bool isSolved() const;
        const std::vector<double>& getEigenValues() const;
        const Matrix<std::complex<double>>& getEigenVectors() const;
};


#endif //HERMITIAN