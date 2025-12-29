#ifndef AURELIA_FOUNDATION_MATRIX_H
#define AURELIA_FOUNDATION_MATRIX_H

#include <array>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <type_traits>
#include "Vector.h" 

namespace Aurelia {
namespace Math {


    template <typename T, size_t Rows, size_t Cols>
    class Matrix {
        static_assert(std::is_floating_point<T>::value, 
            "Matrix<T> requires floating point types.");

    private:
        std::array<T, Rows * Cols> data_;


        constexpr size_t idx(size_t r, size_t c) const {
            return r * Cols + c;
        }

    public:

        Matrix() { data_.fill(T(0)); }

        explicit Matrix(T val) { data_.fill(val); }


        Matrix(std::initializer_list<std::initializer_list<T>> list) {
            size_t r = 0;
            for (const auto& row_list : list) {
                size_t c = 0;
                for (const auto& val : row_list) {
                    if (r < Rows && c < Cols) data_[idx(r, c)] = val;
                    c++;
                }
                r++;
            }
        }

        Matrix(const Matrix&) = default;
        Matrix& operator=(const Matrix&) = default;


        
        size_t rows() const { return Rows; }
        size_t cols() const { return Cols; }

        T& operator()(size_t r, size_t c) {
            #ifdef DEBUG_MATH
            if (r >= Rows || c >= Cols) throw std::out_of_range("Matrix(): Index out of bounds");
            #endif
            return data_[idx(r, c)];
        }

        const T& operator()(size_t r, size_t c) const {
            #ifdef DEBUG_MATH
            if (r >= Rows || c >= Cols) throw std::out_of_range("Matrix(): Index out of bounds");
            #endif
            return data_[idx(r, c)];
        }



        Matrix operator+(const Matrix& other) const {
            Matrix res;
            for (size_t i = 0; i < data_.size(); ++i) res.data_[i] = data_[i] + other.data_[i];
            return res;
        }

        Matrix operator-(const Matrix& other) const {
            Matrix res;
            for (size_t i = 0; i < data_.size(); ++i) res.data_[i] = data_[i] - other.data_[i];
            return res;
        }

        Matrix operator*(T scalar) const {
            Matrix res;
            for (size_t i = 0; i < data_.size(); ++i) res.data_[i] = data_[i] * scalar;
            return res;
        }


        template <size_t OtherCols>
        Matrix<T, Rows, OtherCols> operator*(const Matrix<T, Cols, OtherCols>& other) const {
            Matrix<T, Rows, OtherCols> res; 
            
            for (size_t i = 0; i < Rows; ++i) {
                for (size_t k = 0; k < Cols; ++k) {
                    T temp = (*this)(i, k);
                    if (std::abs(temp) < 1e-15) continue; 
                    for (size_t j = 0; j < OtherCols; ++j) {
                        res(i, j) += temp * other(k, j);
                    }
                }
            }
            return res;
        }


        Vector<T, Rows> operator*(const Vector<T, Cols>& vec) const {
            Vector<T, Rows> res;
            for (size_t i = 0; i < Rows; ++i) {
                T sum = T(0);
                for (size_t j = 0; j < Cols; ++j) {
                    sum += (*this)(i, j) * vec[j];
                }
                res[i] = sum;
            }
            return res;
        }



        Matrix<T, Cols, Rows> transpose() const {
            Matrix<T, Cols, Rows> res;
            for (size_t i = 0; i < Rows; ++i) {
                for (size_t j = 0; j < Cols; ++j) {
                    res(j, i) = (*this)(i, j);
                }
            }
            return res;
        }

        T trace() const {
            static_assert(Rows == Cols, "Trace requires square matrix.");
            T sum = T(0);
            for (size_t i = 0; i < Rows; ++i) sum += (*this)(i, i);
            return sum;
        }

        static Matrix identity() {
            static_assert(Rows == Cols, "Identity requires square matrix.");
            Matrix res;
            for (size_t i = 0; i < Rows; ++i) res(i, i) = T(1);
            return res;
        }

        T determinant() const {
            static_assert(Rows == Cols, "Determinant requires square matrix.");
            
            if constexpr (Rows == 3) {

                return  data_[0] * (data_[4] * data_[8] - data_[5] * data_[7]) -
                        data_[1] * (data_[3] * data_[8] - data_[5] * data_[6]) +
                        data_[2] * (data_[3] * data_[7] - data_[4] * data_[6]);
            } else {

                return T(0); 
            }
        }

        Matrix inverse() const {
            static_assert(Rows == Cols, "Inverse requires square matrix.");
            
            if constexpr (Rows == 3) {

                T det = determinant();
                if (std::abs(det) < 1.0e-12) throw std::runtime_error("Matrix Inverse: Singular.");
                
                T invDet = T(1) / det;
                Matrix res;

                res(0,0) = (data_[4]*data_[8] - data_[5]*data_[7]) * invDet;
                res(0,1) = (data_[2]*data_[7] - data_[1]*data_[8]) * invDet;
                res(0,2) = (data_[1]*data_[5] - data_[2]*data_[4]) * invDet;

                res(1,0) = (data_[5]*data_[6] - data_[3]*data_[8]) * invDet;
                res(1,1) = (data_[0]*data_[8] - data_[2]*data_[6]) * invDet;
                res(1,2) = (data_[2]*data_[3] - data_[0]*data_[5]) * invDet;

                res(2,0) = (data_[3]*data_[7] - data_[4]*data_[6]) * invDet;
                res(2,1) = (data_[1]*data_[6] - data_[0]*data_[7]) * invDet;
                res(2,2) = (data_[0]*data_[4] - data_[1]*data_[3]) * invDet;

                return res;
            } else {
                return Matrix();
            }
        }

        friend std::ostream& operator<<(std::ostream& os, const Matrix& m) {
            os << std::scientific << std::setprecision(4);
            for (size_t i = 0; i < Rows; ++i) {
                os << "[ ";
                for (size_t j = 0; j < Cols; ++j) {
                    os << std::setw(12) << m(i, j) << " ";
                }
                os << "]\n";
            }
            return os;
        }
    };

} 
} 

#endif 