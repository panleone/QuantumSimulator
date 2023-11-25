#include "grid.h"
#include <cassert>

Grid::Grid(double leftBound, double rightBound, std::size_t dim) : leftBound(leftBound), rightBound(rightBound), dim(dim) {
    assert(dim >= 2 && "Grid must contain AT LEAST 2 points");
    assert(rightBound > leftBound && "rightBound must be bigger than leftBound!");
}
std::size_t Grid::getDim() const {
    return dim;
}
double Grid::getLeftBound() const {
    return leftBound;
}
double Grid::getRightBound() const {
    return rightBound;
}
double Grid::getNthPoint(std::size_t n) const {
    assert(n < dim && "the grid has only dim points!");
    return leftBound + n*getStep();
}
double Grid::getStep() const {
    return (rightBound - leftBound)/(dim - 1.0);
}

std::ostream& operator<<(std::ostream& o, const Grid& grid){
    o << "[";
    for(size_t i = 0; i < grid.dim; i++ ){
        o << " " << grid.getNthPoint(i);
        if(i != grid.dim -1){
            o << ",";
        }
    }
    o << "]";
    return o;
}
