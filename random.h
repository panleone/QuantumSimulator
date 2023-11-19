#ifndef RANDOM
#define RANDOM

#include <random>
#include <complex>
#include "matrix.h"
#include "hermitianMatrix.h"
/**
 * @brief Generate a random complex number
 * 
 * @param lowerBound - Minimum value generated as real/imaginary part
 * @param upperBound - Maximum Value generated as real/imaginary part
 * @param isReal - If true the output will be real
 * @return std::complex<double> 
 */
std::complex<double> randomComplexNumber(double lowerBound = -1.0, double upperBound = 1.0, bool isReal = false);
/**
 * @brief Generate a random complex hermitian matrix, only lower triangular values are stored!
 * 
 * @param dim - The dimension of the matrix
 * @return Matrix<std::complex<double>> 
 */
HermitianMatrix randomComplexLTHermitianMatrix(std::size_t dim);
/**
 * @brief Generate a random complex matrix of given dimension
 * 
 * @param dimX - The number of rows
 * @param dimY - The number of columns
 * @return HermitianMatrix
 */
Matrix<std::complex<double>> randomComplexMatrix(std::size_t dimX, std::size_t dimY);

#endif // RANDOM