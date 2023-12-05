#include "./multiParticleWF.h"
#include <cassert>
#include <iostream>

MultiParticleWF::MultiParticleWF(size_t dim, size_t N) : dim(dim), N(N), Matrix<std::complex<double>>(pow(dim,N),1) {};

std::complex<double>& MultiParticleWF::operator()(size_t i) {
    return Matrix<std::complex<double>>::operator()(i,0);
}

const std::complex<double>& MultiParticleWF::operator()(size_t i) const {
    assert(i < pow(dim,N) && "MultiParticleWF element out of bound!");
    return Matrix<std::complex<double>>::operator()(i,0);
}

double MultiParticleWF::getNorm() const {
    double sum = 0;
    for(const auto& element : this->matData){
        sum += std::norm(element);
    }
    return sqrt(sum);
}

void MultiParticleWF::normalize() {
    (*this)/=getNorm();
}

MultiParticleSeparableWF::MultiParticleSeparableWF(size_t dim, size_t N){
    this->N = N;
    this->dim = dim;
    for(int i = 0; i < N ; i++){
        this->wfData.push_back(SingleParticleWF(dim));
    }
};

std::complex<double>& MultiParticleSeparableWF::operator()(size_t i, size_t j){
    assert(i < N && "MultiParticleSeparableWF has only N particles!");
    return this->wfData.at(i)(j);
}

const std::complex<double>& MultiParticleSeparableWF::operator()(size_t i, size_t j) const{
    assert(i < N && "MultiParticleSeparableWF has only N particles!");
    return this->wfData.at(i)(j);
}

void MultiParticleSeparableWF::normalize() {
    for(auto& wf : wfData){
        wf.normalize();
    }
}