#ifndef AURELIA_FOUNDATION_VECTOR_H
#define AURELIA_FOUNDATION_VECTOR_H

#include <vector>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <initializer_list>
#include <type_traits> 
#include <numeric>     


#include "Matrix.h"

namespace Aurelia {
namespace Math {

    template <typename T>
    class Vector {
        static_assert(std::is_floating_point<T>::value, 
            "Vector<T> requires floating point types (float, double, long double).");

    private:
        size_t size_;
        std::vector<T> data_;

    public:
        Vector() : size_(0) {}

        explicit Vector(size_t s) : size_(s), data_(s, T(0)) {}

        Vector(size_t s, T val) : size_(s), data_(s, val) {}

        Vector(std::initializer_list<T> list) : size_(list.size()), data_(list) {}

        Vector(const std::vector<T>& v) : size_(v.size()), data_(v) {}

        Vector(const Vector&) = default;
        Vector(Vector&&) noexcept = default;
        Vector& operator=(const Vector&) = default;
        Vector& operator=(Vector&&) noexcept = default;
    
        using iterator = typename std::vector<T>::iterator;
        using const_iterator = typename std::vector<T>::const_iterator;

        iterator begin() { return data_.begin(); }
        iterator end() { return data_.end(); }
        const_iterator begin() const { return data_.begin(); }
        const_iterator end() const { return data_.end(); }
        const_iterator cbegin() const { return data_.cbegin(); }
        const_iterator cend() const { return data_.cend(); }


        size_t size() const { return size_; }

        T& operator[](size_t i) {
            #ifdef DEBUG_MATH
            if (i >= size_) throw std::out_of_range("Vector[]: Index out of bounds");
            #endif
            return data_[i];
        }

        const T& operator[](size_t i) const {
            #ifdef DEBUG_MATH
            if (i >= size_) throw std::out_of_range("Vector[]: Index out of bounds");
            #endif
            return data_[i];
        }


        Vector operator+(const Vector& other) const {
            if (size_ != other.size_) throw std::invalid_argument("Vector +: Dimension mismatch");
            Vector res(size_);
            for (size_t i = 0; i < size_; ++i) res.data_[i] = data_[i] + other.data_[i];
            return res;
        }

        Vector operator-(const Vector& other) const {
            if (size_ != other.size_) throw std::invalid_argument("Vector -: Dimension mismatch");
            Vector res(size_);
            for (size_t i = 0; i < size_; ++i) res.data_[i] = data_[i] - other.data_[i];
            return res;
        }


        Vector operator*(T scalar) const {
            Vector res(size_);
            for (size_t i = 0; i < size_; ++i) res.data_[i] = data_[i] * scalar;
            return res;
        }


        friend Vector operator*(T scalar, const Vector& v) {
            return v * scalar; 
        }


        Vector operator-() const {
            Vector res(size_);
            for (size_t i = 0; i < size_; ++i) res.data_[i] = -data_[i];
            return res;
        }

        Vector& operator+=(const Vector& other) {
            if (size_ != other.size_) throw std::invalid_argument("Vector +=: Dimension mismatch");
            for (size_t i = 0; i < size_; ++i) data_[i] += other.data_[i];
            return *this;
        }

        Vector& operator-=(const Vector& other) {
            if (size_ != other.size_) throw std::invalid_argument("Vector -=: Dimension mismatch");
            for (size_t i = 0; i < size_; ++i) data_[i] -= other.data_[i];
            return *this;
        }

        Vector& operator*=(T scalar) {
            for (size_t i = 0; i < size_; ++i) data_[i] *= scalar;
            return *this;
        }


        T dot(const Vector& other) const {
            if (size_ != other.size_) throw std::invalid_argument("Vector dot: Dimension mismatch");
            T sum = T(0);
            for (size_t i = 0; i < size_; ++i) sum += data_[i] * other.data_[i];
            return sum;
        }


        Vector cross(const Vector& other) const {
            if (size_ != 3 || other.size_ != 3) 
                throw std::runtime_error("Vector cross: Strictly defined for 3D vectors.");
            

            Vector res(3); 
            res.data_[0] = data_[1] * other.data_[2] - data_[2] * other.data_[1];
            res.data_[1] = data_[2] * other.data_[0] - data_[0] * other.data_[2];
            res.data_[2] = data_[0] * other.data_[1] - data_[1] * other.data_[0];
            return res;
        }

        Matrix<T> outer(const Vector& other) const {
            Matrix<T> res(size_, other.size());
            for (size_t i = 0; i < size_; ++i) {
                for (size_t j = 0; j < other.size(); ++j) {
                    res(i, j) = data_[i] * other[j];
                }
            }
            return res;
        }

 

        T squaredNorm() const {
            T sum = T(0);
            for (const auto& val : data_) sum += val * val;
            return sum;
        }

        T norm() const {
            return std::sqrt(squaredNorm());
        }

        Vector normalized() const {
            T n = norm();
            if (n < Aurelia::Config::NUMERICAL_EPSILON) 
                throw std::runtime_error("Vector normalized: Singularity (Zero vector).");
            return (*this) * (T(1) / n);
        }

        void normalize() {
            T n = norm();
            if (n < Aurelia::Config::NUMERICAL_EPSILON) 
                throw std::runtime_error("Vector normalize: Singularity (Zero vector).");
            T inv = T(1) / n;
            for (auto& val : data_) val *= inv;
        }


        friend std::ostream& operator<<(std::ostream& os, const Vector& v) {
            os << "(";
            for (size_t i = 0; i < v.size_; ++i) {
                os << (i == 0 ? "" : ", ") << v.data_[i];
            }
            os << ")";
            return os;
        }
    };


    using Vector3 = Vector<long double>;
    using VectorX = Vector<long double>;

} 
} 

#endif