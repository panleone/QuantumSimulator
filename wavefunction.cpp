#include "wavefunction.h"
#include "matrix.h"
#include <cassert>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <fftw3.h>

double WaveFunction::getNorm() const {
    double sum = 0;
    for(const auto& element : this->matData){
        sum += std::norm(element);
    }
    return sqrt(sum);
}
//TODO: recheck this constructor
WaveFunction::WaveFunction(const Grid& spaceGrid) : Matrix<std::complex<double>>(spaceGrid.getDim(), 1), spaceGrid(spaceGrid), momentumGrid(spaceGrid) 
{
    size_t dim = spaceGrid.getDim();
    // TODO: I'm not 100% sure if the current way I generate the momentum grid
    // works also for odd values of dim, recheck and in case remove the assert
    assert(dim%2 == 0 && "Wavefunction's grid must have an even number of elements");
    momentumGrid = Grid(-M_PI*(dim)/(spaceGrid.getStep()*dim), M_PI*(dim-2)/(spaceGrid.getStep()*dim), dim);
    fourierTransformed = false;
}

bool WaveFunction::isFourierTransformed() const {
    return fourierTransformed;
}

// TODO: in this way each time we want to FFT a new plan is created
// a better way would be creating the plane once and using it many times
void WaveFunction::FFT() {
    fftw_plan p;
    const size_t n = spaceGrid.getDim();
    // Whether we should transform or anti transform
    const int8_t sign = fourierTransformed ? FFTW_BACKWARD : FFTW_FORWARD;
    p = fftw_plan_dft_1d(n,
                     reinterpret_cast<fftw_complex*>(&*matData.begin()),
                     reinterpret_cast<fftw_complex*>(&*matData.begin()),
                     sign, FFTW_ESTIMATE);
    fftw_execute(p);
    // FFTW does not normalize so we have do divide by the sqrt of the dimension
    // it is a O(N) operation so it's ok
    (*this)/=sqrt(n);
    fourierTransformed = !fourierTransformed;
    fftw_destroy_plan(p);
}

double WaveFunction::getMomentum(size_t i) const {
    const size_t dim = momentumGrid.getDim();
    assert(i < dim && "Trying to access out of bound elements");
    return momentumGrid.getNthPoint((i+dim/2) % dim);
}
double WaveFunction::getPosition(size_t i) const {
    return spaceGrid.getNthPoint(i);
}
std::complex<double>& WaveFunction::operator()(std::size_t x) {
    return Matrix<std::complex<double>>::operator()(x,0);
}

const std::complex<double>& WaveFunction::operator()(std::size_t x) const {
    return Matrix<std::complex<double>>::operator()(x,0);
}

double WaveFunction::positionAverage() const {
    assert(!fourierTransformed);
    double posAverage = 0;
    for(size_t i = 0; i < spaceGrid.getDim(); i++){
        posAverage += spaceGrid.getNthPoint(i)*std::norm((*this)(i));
    }
    return posAverage/pow(getNorm(),2);
}

double WaveFunction::positionVariance() const {
    assert(!fourierTransformed);
    double posAverage = 0;
    double posSquareAverage = 0;
    for(size_t i = 0; i < spaceGrid.getDim(); i++){
        double xi = spaceGrid.getNthPoint(i);
        posAverage += xi*std::norm((*this)(i));
        posSquareAverage += pow(xi,2)*std::norm((*this)(i));
    }
    double  normSquared = pow(getNorm(),2);
    return (posSquareAverage/normSquared - pow(posAverage/normSquared,2));
}

size_t WaveFunction::getDim() const {
    return getDimX();
}