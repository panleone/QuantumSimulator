#ifndef HERMITIAN
#define HERMITIAN

#include "matrix.h"
#include <algorithm>
#include <complex>
#include <cstddef>
#include <vector>

class ComplexVector;
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
        HermitianMatrix(SquareMatrix<std::complex<double>>&& mat) : SquareMatrix<std::complex<double>>(std::move(mat)), eigenValues(getDim()) {};
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
        // Returns the N-th eigenvector of the matrix
        const ComplexVector getNthEigenVector(size_t N) const;
};

class ComplexMatrix : public Matrix<std::complex<double>>
{
    public:
        ComplexMatrix(const Matrix<std::complex<double>>& mat) : Matrix<std::complex<double>>(mat) {};
        ComplexMatrix(const SquareMatrix<std::complex<double>>& mat) : Matrix<std::complex<double>>(mat) {};
        ComplexMatrix(const SquareMatrix<std::complex<double>>&& mat) : Matrix<std::complex<double>>(std::move(mat)) {};
        explicit ComplexMatrix(std::size_t nRows, std::size_t nColumns) : Matrix<std::complex<double>>(nRows,nColumns){};
        ComplexMatrix getAdjoint() const;
};

class ComplexVector : public ComplexMatrix
{
    public:
        explicit ComplexVector(std::size_t dim) : ComplexMatrix(dim,1){};
        std::complex<double>& operator()(size_t i) {return Matrix<std::complex<double>>::operator()(i,0);};
        const std::complex<double>& operator()(size_t i) const {return Matrix<std::complex<double>>::operator()(i,0);};
        
        size_t getDim() const {return getDimX();};
        HermitianMatrix getProjector() const;
    private:
        // Those are private since we use getDim() in place
        using Matrix<std::complex<double>>::getDimX;
        using Matrix<std::complex<double>>::getDimY;
};

#endif //HERMITIAN