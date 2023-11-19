#include "random.h"
#include <chrono>  // for std::chrono::system_clock::now()

std::complex<double> randomComplexNumber(double lowerBound, double upperBound, bool isReal){
   std::uniform_real_distribution<double> unif(lowerBound, upperBound);
   std::default_random_engine re;
   re.seed(std::chrono::system_clock::now().time_since_epoch().count());
   return std::complex<double>(unif(re), isReal ? 0.0 : unif(re));
}

HermitianMatrix randomComplexLTHermitianMatrix(std::size_t dim){
    HermitianMatrix randomM = HermitianMatrix(dim);
    for (size_t i=0; i<dim; i++){
        for(size_t j=i; j<dim; j++){
            auto rComplex = randomComplexNumber(-100.0, 100.0, i==j);
            randomM(j,i) = rComplex;
        }
    }
    return randomM;
}

Matrix<std::complex<double>> randomComplexMatrix(std::size_t dimX, std::size_t dimY){
    Matrix<std::complex<double>> randomM(dimX, dimY);
    for (size_t i=0; i<dimX; i++){
        for(size_t j=0; j<dimY; j++){
            auto rComplex = randomComplexNumber(-100.0, 100.0);
            randomM(i,j) = rComplex;
        }
    }
    return randomM;
}