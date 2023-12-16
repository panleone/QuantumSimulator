#ifndef MATRIX
#define MATRIX

#include <cassert> // assert
#include <cstddef> // std::size_t
#include <iostream>
#include <ostream>
#include <vector>

template<typename T>
class Matrix;
template<typename T>
std::ostream& operator<<(std::ostream& o, const Matrix<T>& m1);
template<typename T>
Matrix<T> operator-(const Matrix<T>& m1, const Matrix<T>& m2);
template<typename T>
Matrix<T> operator+(const Matrix<T>& m1, const Matrix<T>& m2);
template<typename T>
Matrix<T> operator*(const Matrix<T>& m1, const T& c);
template<typename T>
Matrix<T> operator*(const T& c ,const Matrix<T>& m1);
template<typename T>
Matrix<T> tens_product(const Matrix<T>& m1 ,const Matrix<T>& m2);

template<typename T>
/**
 * @brief Template class that represents a generic Matrix M, basic functions are provided
 * 
 */
class Matrix{
    protected:
        // Number of rows in the matrix
        std::size_t nRows;
        // Data stored in the matrix 
        std::vector<T> matData;

    public:
        explicit Matrix(std::size_t nRows, std::size_t nColumns);
        Matrix(const Matrix<T>& m1) = default;

        friend std::ostream& operator<< <> (std::ostream& o,const Matrix<T>& m);
        std::size_t getDimX() const {return nRows;};
        std::size_t getDimY() const {return matData.size()/nRows;};

        T& operator()(std::size_t x, std::size_t y);
        const T& operator()(std::size_t x, std::size_t y) const;

        // Basic operations
        friend Matrix<T> operator+ <>(const Matrix<T>& m1, const Matrix<T>& m2);
        friend Matrix<T> operator- <>(const Matrix<T>& m1, const Matrix<T>& m2);
        friend Matrix<T> operator* <>(const Matrix<T>& m1, const T& c);
        friend Matrix<T> operator* <>(const T& c,const Matrix<T>& m1);
        friend Matrix<T> tens_product <>(const Matrix<T>& m1, const Matrix<T>& m2);

        Matrix<T> operator-() const;
        Matrix<T>& operator+=(const Matrix<T>& m1);
        Matrix<T>& operator-=(const Matrix<T>& m1);
        Matrix<T>& operator*=(const T& c1);
        Matrix<T>& operator/=(const T& c1);
};

template<typename T>
/**
 * @brief Construct a new Matrix< T>:: Matrix object
 * 
 * @param dimX - Number of rows
 * @param dimY - Number of columns 
 */
Matrix<T>::Matrix(std::size_t nRows, std::size_t nColumns):nRows(nRows), matData(nRows*nColumns){}

template<typename T>
/**
 * @brief Overloaded << operator
 * 
 * @param o 
 * @param m1 
 * @return std::ostream& 
 */
std::ostream& operator<<(std::ostream& o, const Matrix<T>& m1){
    //TODO: find a more efficient way to print
    for(size_t j = 0;j < m1.getDimX(); j++){
        for(size_t i = 0; i < m1.getDimY(); i++){       
            o << m1.matData.at(j + i*m1.getDimX()) << " ";
        }
        o << "\n";
    }
    return o;
}

template<typename T>
Matrix<T> operator+(const Matrix<T>& m1, const Matrix<T>& m2){
    assert(m1.getDimX() == m2.getDimX() && m1.getDimY() == m2.getDimY() && "Cannot sum matrices of different sizes!");
    Matrix<T> res = m1;
    return res+=m2;
}

template<typename T>
Matrix<T> Matrix<T>::operator-() const{
    Matrix<T> res = *this;
    for (T& el : res.matData){
        el = -el;
    }
    return res;
}

template<typename T>
Matrix<T> operator-(const Matrix<T>& m1, const Matrix<T>& m2){
    assert(m1.getDimX() == m2.getDimX() && m1.getDimY() == m2.getDimY() && "Cannot subtract matrices of different sizes!");
    return m1+(-m2);
}

template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& m1){
    assert(m1.getDimX() == getDimX() && m1.getDimY() == getDimY() && "Cannot sum matrices of different sizes!");
    for(size_t i = 0; i< getDimX()*getDimY(); i++){
        matData[i] += m1.matData[i];
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& m1){
    assert(m1.getDimX() == getDimX() && m1.getDimY() == getDimY() && "Cannot subtract matrices of different sizes!");
    for(size_t i = 0; i< getDimX()*getDimY(); i++){
        matData[i] -= m1.matData[i];
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(const T& c1){
    for(size_t i = 0; i< getDimX()*getDimY(); i++){
        matData[i] *= c1;
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator/=(const T& c1){
    return (*this)*= (1.0/c1);
}

template<typename T>
/**
 * @brief Get a reference to the element at position (x,y)
 * 
 * @param x - row we want to access
 * @param y - column we want to access
 * @return T& 
 */
T& Matrix<T>::operator()(std::size_t x, std::size_t y) {
    assert( x < getDimX() && "x coordinate out of bound in get()!");
    assert( y < getDimY() && "y coordinate out of bound in get()!");
    return matData.at(x + y*getDimX());
}

template<typename T>
Matrix<T> operator*(const Matrix<T>& m1, const T& c) {
    Matrix<T> res = m1;
    for (T& el : res.matData) {
        el *= c;
    }
    return res;
}

template<typename T>
Matrix<T> operator*(const T& c, const Matrix<T>& m1) {
    return m1*c;
}

template<typename T>
Matrix<T> tens_product(const Matrix<T>& m1 ,const Matrix<T>& m2) {
    Matrix<T> res(m1.getDimX()*m2.getDimX(), m1.getDimY()*m2.getDimY());
    for(int i = 0; i < m1.getDimX(); i++){
        for(int j = 0; j < m1.getDimY(); j++){
            for(int k = 0; k < m2.getDimX(); k++){
                for(int l = 0; l < m2.getDimY(); l++){
                    res(m2.getDimX()*i+k, m2.getDimY()*j+l) = m1(i,j)*m2(k,l);
                }    
            }
        }
    }
    return res;
}

template<typename T>
const T& Matrix<T>::operator()(std::size_t x, std::size_t y) const{
    assert( x < getDimX() && "x coordinate out of bound in get()!");
    assert( y < getDimY() && "y coordinate out of bound in get()!");
    return matData.at(x + y*getDimX());
}

template<typename T>
class SquareMatrix;
template<typename T>
SquareMatrix<T> tens_product(const SquareMatrix<T>& m1 ,const SquareMatrix<T>& m2);

template<typename T>
/**
 * @brief Template class that represents a generic Square Matrix M, basic functions are provided
 * 
 */
class SquareMatrix : public Matrix<T>{
    public:
        explicit SquareMatrix(std::size_t size) : Matrix<T>(size, size){};
        SquareMatrix(const Matrix<T>&& m) : Matrix<T>(m) {
            assert(m.getDimX() == m.getDimY() && "Matrix is not a square amtrix");
        }
        size_t getDim() const { return getDimX();}
        // Returns the trace of the matrix
        T trace() const;
        friend SquareMatrix<T> tens_product <>(const SquareMatrix<T>& m1, const SquareMatrix<T>& m2);
    private:
        // Those are private since we use getDim() in place
        using Matrix<T>::getDimX;
        using Matrix<T>::getDimY;
};

template<typename T>
T SquareMatrix<T>::trace() const{
    T trace;
    for(std::size_t i = 0; i < getDim(); i++){
        trace += (*this)(i,i);
    }
    return trace;
}

template<typename T>
SquareMatrix<T> tens_product(const SquareMatrix<T>& m1 ,const SquareMatrix<T>& m2) {
    return tens_product(static_cast<const Matrix<T>&>(m1), static_cast<const Matrix<T>&>(m2));
}

#endif // MATRIX
        

