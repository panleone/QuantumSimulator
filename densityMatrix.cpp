#include "./densityMatrix.h"
#include "hermitianMatrix.h"
#include "matrix.h"
#include <cassert>
#include <cmath>
#include <cstddef>

PureDensityMatrix::PureDensityMatrix(const MultiParticleWF& waveFunction) : HermitianMatrix(waveFunction.getDim()){
    double sq_norm = pow(waveFunction.getNorm(), 2);
    N = waveFunction.getN();
    for(size_t i = 0; i < waveFunction.getDim(); i++){
        for(size_t j = 0; j < waveFunction.getDim(); j++){
           (*this)(i,j) = std::conj(waveFunction(i))*waveFunction(j)/sq_norm;
        }
    }
}

DensityMatrix::DensityMatrix(size_t size, size_t N) : HermitianMatrix(size), N(N){};

DensityMatrix::DensityMatrix(const std::vector<PureDensityMatrix>& pureDensityMatrices, const std::vector<double>& probabilities) : HermitianMatrix(pureDensityMatrices.at(0).getDim()){
    size_t dimension = pureDensityMatrices.at(0).getDim();
    N = pureDensityMatrices.at(0).getN();
    // Sanity check on the input
    assert(pureDensityMatrices.size() == probabilities.size() && "DensityMatrix: size of vectors does not match!");
    for(const auto& pureDensityMatrix: pureDensityMatrices){
        assert(pureDensityMatrix.getDim() == dimension && "DensityMatrix: dimension of pureDensityMatrices does not match!");
        assert(pureDensityMatrix.getN() == N && "DensityMatrix: number of particles of pureDensityMatrices does not match!");
    }

    // Fill the matrix
    for(size_t i = 0; i < dimension; i++){
        for(size_t j = 0; j < dimension; j++){
            for(size_t k = 0; k < pureDensityMatrices.size(); k++){
                const auto& pureDensityMatrix = pureDensityMatrices.at(k);
                const auto& probability = probabilities.at(k);
                (*this)(i,j) += pureDensityMatrix(i,j)*probability;
            }
        }
    }
}

DensityMatrix::DensityMatrix(const std::vector<MultiParticleWF>& waveFunctions, const std::vector<double>& probabilities) : HermitianMatrix(waveFunctions.at(0).getDim()){
    size_t dimension = waveFunctions.at(0).getDim();
    N = waveFunctions.at(0).getN();
    // Sanity check on the input
    assert(waveFunctions.size() == probabilities.size() && "DensityMatrix: size of vectors does not match!");
    for(const auto& waveFunction: waveFunctions){
        assert(waveFunction.getDim() == dimension && "DensityMatrix: dimension of wave functions does not match!");
        assert(waveFunction.getN() == N && "DensityMatrix: number of particles of wave functions does not match!");
    }

    std::vector<double> sq_norms(dimension);
    for(size_t k = 0; k < waveFunctions.size(); k++){
        sq_norms.at(k) = pow(waveFunctions.at(k).getNorm(), 2);
    } 
    // Fill the matrix
    for(size_t i = 0; i < dimension; i++){
        for(size_t j = 0; j < dimension; j++){
            for(size_t k = 0; k < waveFunctions.size(); k++){
                const auto& waveFunction = waveFunctions.at(k);
                const auto& probability = probabilities.at(k);
                (*this)(i,j) += std::conj(waveFunction(i))*waveFunction(j)*probability/sq_norms.at(k);
            }
        }
    }
}

DensityMatrix::DensityMatrix(const PureDensityMatrix& pureDensityMatrix) : HermitianMatrix(pureDensityMatrix.getDim()){
    N = pureDensityMatrix.getN();
    for(size_t i = 0; i < getDim(); i++){
        for(size_t j = 0; j < getDim(); j++){ 
            (*this)(i,j) = pureDensityMatrix(i,j);
        }
    }
}

DensityMatrix DensityMatrix::partialTrace(size_t subSystem){
    // At the moment only the case N = 2 is supported
    assert(N == 2 && "At the moment partial trace can only be calculated for a 2 particle system");
    assert(subSystem < 2 && "partialTrace: The system has only 2 particles!");
    // The original DensityMatrix had dimension D^2 x D^2, the new one will have dimension D x D
    size_t redSize = sqrt(getDim());
    DensityMatrix res(redSize, 1);
    for( size_t i = 0; i < redSize; i++){
        for( size_t j = 0; j < redSize; j++){
            for( size_t k = 0; k < redSize; k++){
                if (subSystem == 0){
                    res(i,j) += (*this)(i*redSize + k, j*redSize + k);
                }else{
                    res(i,j) += (*this)(i + k*redSize, j + k*redSize);
                }
            }
        }
    }
    return res;
}