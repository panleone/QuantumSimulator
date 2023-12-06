#ifndef MULTI_PARTICLE_WF
#define MULTI_PARTICLE_WF

#include "matrix.h"
#include <cassert>
#include <complex.h>
#include <cstddef>
#include <vector>

/**
 * @brief Describes a system made by many particles
 */
class MultiParticleWF : public Matrix<std::complex<double>> {
    public:
        /**
         * @brief Construct a new Multi Particle WF object
         * 
         * @param N - the number of particles
         * @param dim - hilbert space dimension of each particle
         */
        explicit MultiParticleWF(size_t dim, size_t N);
        // Returns a reference to the i-th element of the WF
        std::complex<double>& operator()(size_t i);
        const std::complex<double>& operator()(size_t i) const;
        // Returns the dimension of the hilbert space of the MultiParticle wave function
        size_t getDim() const { return pow(dim, N); }
        // Returns the norm of the wave function
        double getNorm() const;
        // Returns the total number of particles described by the W.F.
        size_t getN() const {return N;}
        // Normalize the wave function
        void normalize();

    private:
        // Those are private since we use getDim() in place
        using Matrix<std::complex<double>>::getDimX;
        using Matrix<std::complex<double>>::getDimY;
        // Dimension of the space in which each particle lives
        size_t dim;
        // Number of particles of the system
        size_t N;
};

/**
 * @brief Util class to describe the particular case of single particle wave functions
 */
class SingleParticleWF : public MultiParticleWF{
    public:
        explicit SingleParticleWF(size_t dim) : MultiParticleWF(dim, 1) {};
};


/**
 * @brief Describes a system made by many separable particles
 */
class MultiParticleSeparableWF{
    public:
        /**
         * @brief Construct a new Multi Particle Separable WF object
         * 
         * @param N - the number of particles
         * @param dim - hilbert space dimension of each particle
         */
        explicit MultiParticleSeparableWF(size_t dim, size_t N);
        /**
         * @brief Returns a reference to the j-th component of the i-th particle of the WF
         * 
         * @param i - particle we want to access, 0 <= i < N
         * @param j - component of the particle we want to access, 0 <= j < dim
         * @return std::complex<double>& the j-th component of the i-th particle
         */
        std::complex<double>& operator()(size_t i, size_t j);
        const std::complex<double>& operator()(size_t i, size_t j) const;
        // Normalize the wave function
        void normalize();
    private:
        // Data stored in the WF, this has length dim*N
        std::vector<SingleParticleWF> wfData;
        // Dimension of the space in which each particle lives
        size_t dim;
        // Total number of particles
        size_t N;
};

#endif // MULTI_PARTICLE_WF