#ifndef AURELIA_FOUNDATION_VECTOR_H
#define AURELIA_FOUNDATION_VECTOR_H

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <initializer_list>
#include <type_traits>

namespace Aurelia {
namespace Math {

    template <typename T, size_t N>
    class Vector {
        static_assert(std::is_floating_point<T>::value, 
            "Vector<T, N> requires floating point types.");

    private:

        std::array<T, N> data_;

    public:


        Vector() { data_.fill(T(0)); }


        explicit Vector(T val) { data_.fill(val); }


        Vector(std::initializer_list<T> list) {
            size_t i = 0;
            for (auto val : list) {
                if (i < N) data_[i++] = val;
            }

        }


        Vector(const Vector&) = default;
        Vector& operator=(const Vector&) = default;


        T& operator[](size_t i) {
            #ifdef DEBUG_MATH
            if (i >= N) throw std::out_of_range("Vector[]: Index out of bounds");
            #endif
            return data_[i];
        }

        const T& operator[](size_t i) const {
            #ifdef DEBUG_MATH
            if (i >= N) throw std::out_of_range("Vector[]: Index out of bounds");
            #endif
            return data_[i];
        }

        constexpr size_t size() const { return N; }




        Vector operator+(const Vector& other) const {
            Vector res;
            for (size_t i = 0; i < N; ++i) res.data_[i] = data_[i] + other.data_[i];
            return res;
        }

        Vector operator-(const Vector& other) const {
            Vector res;
            for (size_t i = 0; i < N; ++i) res.data_[i] = data_[i] - other.data_[i];
            return res;
        }

        Vector operator-() const {
            Vector res;
            for (size_t i = 0; i < N; ++i) res.data_[i] = -data_[i];
            return res;
        }

        Vector operator*(T scalar) const {
            Vector res;
            for (size_t i = 0; i < N; ++i) res.data_[i] = data_[i] * scalar;
            return res;
        }


        friend Vector operator*(T scalar, const Vector& v) {
            return v * scalar; 
        }

        Vector& operator+=(const Vector& other) {
            for (size_t i = 0; i < N; ++i) data_[i] += other.data_[i];
            return *this;
        }

        Vector& operator-=(const Vector& other) {
            for (size_t i = 0; i < N; ++i) data_[i] -= other.data_[i];
            return *this;
        }

        Vector& operator*=(T scalar) {
            for (size_t i = 0; i < N; ++i) data_[i] *= scalar;
            return *this;
        }




        T dot(const Vector& other) const {
            T sum = T(0);
            for (size_t i = 0; i < N; ++i) sum += data_[i] * other.data_[i];
            return sum;
        }


        Vector cross(const Vector& other) const {
            static_assert(N == 3, "Cross product is only defined for 3D vectors.");
            Vector res;
            res.data_[0] = data_[1] * other.data_[2] - data_[2] * other.data_[1];
            res.data_[1] = data_[2] * other.data_[0] - data_[0] * other.data_[2];
            res.data_[2] = data_[0] * other.data_[1] - data_[1] * other.data_[0];
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

        void normalize() {
            T n = norm();

            T inv = T(1) / n;
            for (auto& val : data_) val *= inv;
        }

        Vector normalized() const {
            Vector v = *this;
            v.normalize();
            return v;
        }



        friend std::ostream& operator<<(std::ostream& os, const Vector& v) {
            os << "(";
            for (size_t i = 0; i < N; ++i) {
                os << (i == 0 ? "" : ", ") << v.data_[i];
            }
            os << ")";
            return os;
        }
    };

}
} 

#endif 