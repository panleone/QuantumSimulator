#ifndef HAMILTONIAN
#define HAMILTONIAN
#include <cstddef>
#include <functional>
#include <optional>
#include <vector>
#include "hermitianMatrix.h"
#include "grid.h"
#include "matrix.h"
#include "wavefunction.h"
/**
 * @brief This class represents a QM 1D single particle hamiltonian
 */
class Hamiltonian : public HermitianMatrix{
    private:
        // The space grid on which we are solving the shrodinger equation
        Grid spaceGrid;
        // Used only for time dependent hamiltonians
        std::optional<Grid> timeGrid = std::nullopt; 
        size_t currentTimeIndex = 0;

        // Mass of the particle
        const double mass;

        // Function representing the potential (x, mass, time) ->  V(x)
        std::function<double(double, double, double)> potential;

        // Generates the kinetic term of the hamiltonian up to fourth order in the grid step
        void setKinetic();
        // Generates the potential term
        void setPotential();
        // Reset the solved flags, so we can re-diagonalize (useful for a time dependent H)
        void reset();
        /**
         * @brief set the current time of a time dependent hamiltonian
         * 
         * @param newTimeIndex - time index smaller than the dimension of the timeGrid
         * @return true iff the operation was succesful
         */
        bool setTime(size_t newTimeIndex);
        
        // Make this operator private since the Hamiltonian should not be changed externally
        using HermitianMatrix::operator();
    public:
        /**
         * @brief Creates a time independent hamiltonian with given potential and default kynetic term
         * 
         * @param spaceGrid - the space grid on which the Shrodinger equation is solved
         * @param mass - mass of the particle
         * @param potFunction - potential energy V: (x, mass, time) -> V(x, m, time) 
         */
        explicit Hamiltonian(const Grid& spaceGrid, double mass, std::function<double(double, double, double)> potFunction);
        // Creates a time dependent hamiltonian, a timeGrid must be provided
        explicit Hamiltonian(const Grid& spaceGrid, const Grid& timeGrid, double mass, std::function<double(double, double, double)> potFunction);

        /**
         * @brief increment the current hamiltonian time by one unit
         * 
         * @return true iff the operation was succesful 
         * (i.e. H is time dependent and time is not already at it's max value)
         */
        bool incrementTime();

        // Applies the time evolution operator on a given W.F.
        void applyTimeEvolutionOperator(WaveFunction& wavefunction) const;
        // Calculates the average energy <\psi| H |\psi>
        double calculateEnergy(WaveFunction& wavefunction) const;
        /**
         * @brief Find eigenvalues and eigenvectors of the hamiltonian
         * evaluated at the current hamiltonian's time
         * @param JOBZ - char: if set to 'N' will compute eigenvalue only, 
         * if set to 'V' will compute also eigenvectors  
         * @return true iff the operation was successful 
         */
        bool solve(char JOBZ);
};

#endif // HAMILTONIAN