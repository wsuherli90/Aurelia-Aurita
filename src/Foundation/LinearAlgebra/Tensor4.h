#ifndef AURELIA_FOUNDATION_TENSOR4_H
#define AURELIA_FOUNDATION_TENSOR4_H

#include <array>      
#include <algorithm>  
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <numeric>
#include "Matrix.h"

namespace Aurelia {
namespace Math {


    template <typename T, size_t N>
    class Tensor4 {
    private:
        std::array<T, N * N * N * N> data_;
        constexpr inline size_t idx(size_t i, size_t j, size_t k, size_t l) const {
            return ((i * N + j) * N + k) * N + l;
        }

    public:


        Tensor4() { data_.fill(T(0)); }

        explicit Tensor4(T val) { data_.fill(val); }



        T& operator()(size_t i, size_t j, size_t k, size_t l) {
            #ifdef DEBUG_MATH
            if (i >= N || j >= N || k >= N || l >= N) 
                throw std::out_of_range("Tensor4 bounds check failed");
            #endif
            return data_[idx(i, j, k, l)];
        }

        const T& operator()(size_t i, size_t j, size_t k, size_t l) const {
            #ifdef DEBUG_MATH
            if (i >= N || j >= N || k >= N || l >= N) 
                throw std::out_of_range("Tensor4 bounds check failed");
            #endif
            return data_[idx(i, j, k, l)];
        }

        void setZero() { data_.fill(T(0)); }
        Matrix<T> doubleContract(const Matrix<T>& F) const {
            if (F.rows() != N || F.cols() != N)
                throw std::invalid_argument("Tensor4: Matrix dimension mismatch.");

            Matrix<T> H(N, N, T(0));
            for (size_t a = 0; a < N; ++a) {
                for (size_t b = 0; b < N; ++b) {
                    T sum = T(0);
                    for (size_t c = 0; c < N; ++c) {
                        for (size_t d = 0; d < N; ++d) {
                            sum += (*this)(a, b, c, d) * F(c, d);
                        }
                    }
                    H(a, b) = sum * T(0.5); 
                }
            }
            return H;
        }

        static Tensor4<T, N> LeviCivita() {
            Tensor4<T, N> eps;
            if constexpr (N == 4) { 
                std::array<int, 4> p = {0, 1, 2, 3};
                do {
                    int inversions = 0;
                    for (size_t i = 0; i < 4; ++i)
                        for (size_t j = i + 1; j < 4; ++j)
                            if (p[i] > p[j]) inversions++;
                    
                    T sign = (inversions % 2 == 0) ? T(1) : T(-1);
                    eps(p[0], p[1], p[2], p[3]) = sign;

                } while (std::next_permutation(p.begin(), p.end()));
            }
            return eps;
        }

        void enforceMajorSymmetry() {
            for (size_t i = 0; i < N; ++i) {
                for (size_t j = 0; j < N; ++j) {
                    for (size_t k = 0; k < N; ++k) {
                        for (size_t l = 0; l < N; ++l) {
                            size_t id1 = idx(i, j, k, l);
                            size_t id2 = idx(k, l, i, j);
                            if (id1 < id2) {
                                T avg = (data_[id1] + data_[id2]) * T(0.5);
                                data_[id1] = data_[id2] = avg;
                            }
                        }
                    }
                }
            }
        }


        void enforceElectromagneticSymmetry() {
            for (size_t c = 0; c < N; ++c)
                for (size_t d = 0; d < N; ++d)
                    for (size_t a = 0; a < N; ++a)
                        for (size_t b = a + 1; b < N; ++b) {
                            T val = ((*this)(a, b, c, d) - (*this)(b, a, c, d)) * T(0.5);
                            (*this)(a, b, c, d) = val;
                            (*this)(b, a, c, d) = -val;
                        }

            for (size_t a = 0; a < N; ++a)
                for (size_t b = 0; b < N; ++b)
                    for (size_t c = 0; c < N; ++c)
                        for (size_t d = c + 1; d < N; ++d) {
                            T val = ((*this)(a, b, c, d) - (*this)(a, b, d, c)) * T(0.5);
                            (*this)(a, b, c, d) = val;
                            (*this)(a, b, d, c) = -val;
                        }
        }


        Tensor4 operator+(const Tensor4& other) const {
            Tensor4 res;
            for (size_t i = 0; i < data_.size(); ++i) {
                res.data_[i] = data_[i] + other.data_[i];
            }
            return res;
        }

        Tensor4 operator*(T scalar) const {
            Tensor4 res;
            for (size_t i = 0; i < data_.size(); ++i) {
                res.data_[i] = data_[i] * scalar;
            }
            return res;
        }

        friend Tensor4 operator*(T scalar, const Tensor4& t) {
            return t * scalar;
        }

        void print() const {
            for (size_t i = 0; i < N; ++i) {
                for (size_t j = 0; j < N; ++j) {
                    std::cout << "Slice (" << i << "," << j << ",:,:):\n";
                    for (size_t k = 0; k < N; ++k) {
                        std::cout << "[ ";
                        for (size_t l = 0; l < N; ++l) {
                            std::cout << std::setw(8) << std::setprecision(3) 
                                      << (*this)(i, j, k, l) << " ";
                        }
                        std::cout << "]\n";
                    }
                    std::cout << "\n";
                }
            }
        }
    };


    using ElasticityTensor = Tensor4<long double, 3>;
    using ConstitutiveTensor = Tensor4<long double, 4>;

} 
} 

#endif 