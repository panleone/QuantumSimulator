#ifndef WAVEFUNCTION 
#define WAVEFUNCTION

#include <complex>
#include <cstddef>
#include <fftw3.h>

#include "matrix.h"
#include "grid.h"
/**
 * @brief This class represents a WaveFunction
 */
class WaveFunction : public Matrix<std::complex<double>> {
    public:
        /**
         * @brief Construct a new Wave Function 
         * 
         * @param spaceGrid - This is the spatial grid where the wavefunction lives, 
         * at the moment only grids with an even number of points are accepted
         */
        WaveFunction(const Grid& spaceGrid);
        // Returns the norm of the W.F.
        double getNorm() const;

        //return true iff the W.F. is represented in momentum space
        bool isFourierTransformed() const;
        /**
         * @brief Performs fast fourier transform or anti-transform 
         * depending on the current state of the wave function
         */
        void FFT();

        // Get momentum corresponding to the i-th component of the W.F.
        double getMomentum(size_t i) const;
        // Get position corresponding to the i-th component of the W.F.
        double getPosition(size_t i) const;

        // Returns the average position of the W.F.
        double positionAverage() const;
        // Returns the variance on the average position of the W.F.
        double positionVariance() const;
        // Returns the dimension of the vector representing the W.F.
        size_t getDim() const;

        //Operator i -> i-th component of the W.F.
        std::complex<double>& operator()(std::size_t x);
        const std::complex<double>& operator()(std::size_t x) const;
    private:
        // Those are private since we use getDim() in place
        using Matrix<std::complex<double>>::getDimX;
        using Matrix<std::complex<double>>::getDimY;
        // Grid on which the wave function is defined
        Grid spaceGrid;
        // Corresponding fourier transformed grid
        Grid momentumGrid;
        // is the wave function represented in momentum space?
        bool fourierTransformed;        
};
#endif //WAVEFUNCTION