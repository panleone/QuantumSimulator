#ifndef HAMILTONIAN
#define HAMILTONIAN
#include <vector>
#include "hermitianMatrix.h"
#include "grid.h"
#include "matrix.h"
/**
 * @brief This class represents a QM 1D single particle hamiltonian
 */

class Hamiltonian{
    private:

        // The grid on which we are solving the shrodinger equation
        Grid grid;

        // Hermitian matrix storing the total hamiltonian
        HermitianMatrix H;

        // Mass of the particle
        const double mass;

        // Generates the kynetic term of the hamiltonian up to second order in the grid step
        void SetKynetic2Ord();
        // Generate the kynetic term of the hamiltonian up to fourth order in the grid step
        void SetKynetic4Ord();
    public:
        /**
         * @brief Creates an hamiltonian with empty potential and default kynetic term
         * 
         * @param dim - The interval [leftBound, rightBound] is discretized in dim points. Must be >= 2.
         * @param leftBound - left extreme of the interval in which we want to solve Shrodinger eq.
         * @param rightBound - right extreme of the interval in which we want to solve Shrodinger eq.
         * @param mass - mass of the particle
         */
        explicit Hamiltonian(std::size_t dim, double leftBound, double rightBound, double mass);
        explicit Hamiltonian(const Grid& grid, double mass);
        /**
         * @brief Set the Potential Pot given the potential function
         * 
         * @param potFunction - function : (x, mass) -> V(x,m) 
         */
        void SetPotential(double (*potFunction)(double, double));
        const HermitianMatrix& getHamiltonian() const;

        /**
         * @brief Find eigenvalues and eigenvectors of the hamiltonian:
         * AFTER CALLING THIS FUNCTION THE HAMILTONIAN IS "DESTROYED"
         * in particular if eigenvalues are computed it's content is replaced
         * with the eigenvalues matrix
         * @param JOBZ - char: if set to 'N' will compute eigenvalue only, 
         * if set to 'V' will compute also eigenvectors  
         * @return true iff the operation was successful 
         */
        bool solve(char JOBZ);
};

#endif // HAMILTONIAN