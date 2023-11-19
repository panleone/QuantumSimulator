#ifndef GRID
#define GRID

#include <cstddef> 

/**
 * @brief This class represents a Grid that discretizes a given real interval [a, b]
 */
class Grid{
    private:
        double leftBound;
        double rightBound;
        // Number of points in the grid
        std::size_t dim;
    public:
        /**
         * @brief Construct the grid [left_bound, right_bound]
         * divided in dim - 1 equally spaced intervals
         * 
         * @param leftBound - left interval
         * @param rightBound - right interval
         * @param dim - number of points >= 2 in the grid
         */
        explicit Grid(double leftBound, double rightBound, std::size_t dim);
        // Getters
        std::size_t getDim() const;
        double getLeftBound() const;
        double getRightBound() const;
        /**
         * @brief Returns the n-th point of the grid
         * 
         * @param n - must be in the range [0, dim)
         * @return double - the coordinate of the n-th point
         */
        double getNthPoint(std::size_t n) const;
        /**
         * @brief Returns the step of the grid
         */
        double getStep() const;
};

#endif //GRID