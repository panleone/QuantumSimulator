#ifndef DENSITY_MATRIX
#define DENSITY_MATRIX

#include <complex>
#include <cstddef>

#include "./hermitianMatrix.h"
#include "multiParticleWF.h"

/**
 * @brief Describes a pure density matrix of any system
 */
class PureDensityMatrix : public HermitianMatrix{
    private:
        // Number of particles described by the Pure Density Matrix
        size_t N;
    public:
        /**
         * @brief Construct a PureDensityMatrix from a waveFanction
         */
        explicit PureDensityMatrix(const MultiParticleWF& waveFunction);
        // Returns the total number of particles described by the density matrix.
        size_t getN() const {return N;}
};

/**
 * @brief Describes a general density matrix
 */
class DensityMatrix : public HermitianMatrix{
    private:
        // Number of particles described by the Density Matrix
        size_t N;
    public:
        /**
         * @brief Construct a DensityMatrix from many PureDensityMatrices
         * 
         * @param pureDensityMatrices - vector of pure density matrices
         * @param probabilities - vector of probabilities, they must sum to 1
         */
        explicit DensityMatrix(const std::vector<PureDensityMatrix>& pureDensityMatrices, const std::vector<double>& probabilities);
        /**
         * @brief Construct a DensityMatrix from wavefunctions
         * 
         * @param waveFunctions - vector of wavefunctions
         * @param probabilities - vector of probabilities, they must sum to 1
         */
        explicit DensityMatrix(const std::vector<MultiParticleWF>& waveFunctions, const std::vector<double>& probabilities);
        /**
         * @brief Construct an empty DensityMatrix from a given N and dim
         * 
         * @param size - size of the density matrix
         * @param N - total number of particles
         */
        explicit DensityMatrix(size_t size, size_t N);
        DensityMatrix(const PureDensityMatrix& densityMatrix);
        /**
         * @brief Computes the reduced density matrix for a given subsystem
         * AT THE MOMENT ONLY THE CASE FOR TWO PARTICLES IS SUPPORTED
         * TODO: generalize to N particle trace
         * @param i - must be between 0 <= subSystem < N
         * @return returns the reduced density matrix for the subSystem-th particle
         */
        DensityMatrix partialTrace(size_t subSystem);
        // Returns the total number of particles described by the density matrix.
        size_t getN() const {return N;}
};

#endif // DENSITY_MATRIX