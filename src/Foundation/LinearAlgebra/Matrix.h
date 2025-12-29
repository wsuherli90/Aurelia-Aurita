#ifndef AURELIA_FOUNDATION_MATRIX_H
#define AURELIA_FOUNDATION_MATRIX_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <algorithm> 
#include <type_traits>


#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace Math {

    template <typename T>
    class Matrix {
        // Removed static_assert to allow Complex<T> types for electrodynamics

    private:
        size_t rows_;
        size_t cols_;
        std::vector<T> data_; 
        inline size_t idx(size_t r, size_t c) const {
            return r * cols_ + c;
        }

    public:
        Matrix() : rows_(0), cols_(0) {}
        Matrix(size_t r, size_t c) : rows_(r), cols_(c) {
             data_.resize(r * c, T(0)); 
        }
        Matrix(size_t r, size_t c, const T& val) : rows_(r), cols_(c), data_(r * c, val) {}
        Matrix(std::initializer_list<std::initializer_list<T>> list) {
            rows_ = list.size();
            cols_ = (rows_ > 0) ? list.begin()->size() : 0;
            data_.reserve(rows_ * cols_); 

            for (const auto& row : list) {
                if (row.size() != cols_) 
                    throw std::invalid_argument("Matrix: Jagged initializer list.");
                data_.insert(data_.end(), row.begin(), row.end());
            }
        }

        Matrix(const std::vector<std::vector<T>>& grid) {
            rows_ = grid.size();
            cols_ = (rows_ > 0) ? grid[0].size() : 0;
            data_.reserve(rows_ * cols_);

            for (const auto& row : grid) {
                if (row.size() != cols_) 
                    throw std::invalid_argument("Matrix: Jagged initialization vector.");
                data_.insert(data_.end(), row.begin(), row.end());
            }
        }

        size_t rows() const { return rows_; }
        size_t cols() const { return cols_; }


        T& operator()(size_t r, size_t c) {
            #ifdef DEBUG_MATH
            if (r >= rows_ || c >= cols_) throw std::out_of_range("Matrix index out of bounds");
            #endif
            return data_[idx(r, c)];
        }


        const T& operator()(size_t r, size_t c) const {
            #ifdef DEBUG_MATH
            if (r >= rows_ || c >= cols_) throw std::out_of_range("Matrix index out of bounds");
            #endif
            return data_[idx(r, c)];
        }


        Matrix operator+(const Matrix& other) const {
            if (rows_ != other.rows_ || cols_ != other.cols_)
                throw std::invalid_argument("Matrix +: Dimension mismatch.");
            Matrix res(rows_, cols_);
            for (size_t i = 0; i < data_.size(); ++i) res.data_[i] = data_[i] + other.data_[i];
            return res;
        }


        Matrix operator-(const Matrix& other) const {
            if (rows_ != other.rows_ || cols_ != other.cols_)
                throw std::invalid_argument("Matrix -: Dimension mismatch.");
            Matrix res(rows_, cols_);
            for (size_t i = 0; i < data_.size(); ++i) res.data_[i] = data_[i] - other.data_[i];
            return res;
        }


        Matrix operator*(T scalar) const {
            Matrix res(rows_, cols_);
            for (size_t i = 0; i < data_.size(); ++i) res.data_[i] = data_[i] * scalar;
            return res;
        }


        Matrix operator*(const Matrix& other) const {
            if (cols_ != other.rows_)
                throw std::invalid_argument("Matrix *: Dimension mismatch (Cols vs Rows).");
            
            Matrix res(rows_, other.cols_, T(0));
            for (size_t i = 0; i < rows_; ++i) {
                for (size_t k = 0; k < cols_; ++k) {
                    T temp = (*this)(i, k);
                    for (size_t j = 0; j < other.cols_; ++j) {
                        res(i, j) += temp * other(k, j);
                    }
                }
            }
            return res;
        }


        Matrix transpose() const {
            Matrix res(cols_, rows_);
            for (size_t i = 0; i < rows_; ++i) {
                for (size_t j = 0; j < cols_; ++j) {
                    res(j, i) = (*this)(i, j);
                }
            }
            return res;
        }


        T trace() const {
            if (rows_ != cols_) throw std::runtime_error("Trace: Matrix must be square.");
            T sum = T(0);
            for (size_t i = 0; i < rows_; ++i) sum += (*this)(i, i);
            return sum;
        }


        static Matrix identity(size_t n) {
            Matrix res(n, n, T(0));
            for (size_t i = 0; i < n; ++i) res(i, i) = T(1);
            return res;
        }

        T determinant() const {
            if (rows_ != cols_) throw std::runtime_error("Determinant: Matrix must be square.");
            if (rows_ == 0) return T(0);
            
            Matrix temp = *this; 
            T det = T(1);
            size_t n = rows_;

            for (size_t i = 0; i < n; ++i) {
                size_t pivot = i;
                for (size_t j = i + 1; j < n; ++j) {
                    if (std::abs(temp(j, i)) > std::abs(temp(pivot, i))) pivot = j;
                }

                if (std::abs(temp(pivot, i)) < Aurelia::Config::NUMERICAL_EPSILON) return T(0);


                if (i != pivot) {
                    for (size_t k = 0; k < n; ++k) std::swap(temp(i, k), temp(pivot, k));
                    det *= -1;
                }

                det *= temp(i, i);


                for (size_t j = i + 1; j < n; ++j) {
                    T factor = temp(j, i) / temp(i, i);
                    for (size_t k = i; k < n; ++k) {
                        temp(j, k) -= factor * temp(i, k);
                    }
                }
            }
            return det;
        }


        Matrix inverse() const {
            if (rows_ != cols_) throw std::runtime_error("Inverse: Matrix must be square.");
            
            size_t n = rows_;
            Matrix augmented(n, 2 * n, T(0));

   
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) augmented(i, j) = (*this)(i, j);
                augmented(i, i + n) = T(1);
            }

      
            for (size_t i = 0; i < n; ++i) {
 
                size_t pivot = i;
                for (size_t j = i + 1; j < n; ++j) {
                    if (std::abs(augmented(j, i)) > std::abs(augmented(pivot, i))) pivot = j;
                }

                if (std::abs(augmented(pivot, i)) < Aurelia::Config::NUMERICAL_EPSILON)
                    throw std::runtime_error("Inverse: Matrix is singular (Metric Degeneration).");

                if (i != pivot) {
                    for (size_t k = 0; k < 2 * n; ++k) 
                        std::swap(augmented(i, k), augmented(pivot, k));
                }


                T div = augmented(i, i);
                T invDiv = T(1) / div; 
                for (size_t k = i; k < 2 * n; ++k) augmented(i, k) *= invDiv;


                for (size_t j = 0; j < n; ++j) {
                    if (i != j) {
                        T factor = augmented(j, i);
                        for (size_t k = i; k < 2 * n; ++k) 
                            augmented(j, k) -= factor * augmented(i, k);
                    }
                }
            }


            Matrix res(n, n);
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    res(i, j) = augmented(i, j + n);
                }
            }
            return res;
        }

        friend std::ostream& operator<<(std::ostream& os, const Matrix& m) {
            os << std::scientific << std::setprecision(6);
            for (size_t i = 0; i < m.rows_; ++i) {
                os << "[ ";
                for (size_t j = 0; j < m.cols_; ++j) {
                    os << std::setw(14) << m(i, j) << " ";
                }
                os << "]\n";
            }
            os.unsetf(std::ios_base::floatfield);
            return os;
        }
    };


    using Matrix3 = Matrix<long double>; 
    using Matrix4 = Matrix<long double>; 
    using Matrix6 = Matrix<long double>; 

} 
} 

#endif 