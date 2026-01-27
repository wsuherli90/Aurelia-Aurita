#ifndef AURELIA_FOUNDATION_COMPLEXMATH_H
#define AURELIA_FOUNDATION_COMPLEXMATH_H

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <limits>

#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace Math {

    template <typename T>
    class Complex {
        static_assert(std::is_floating_point<T>::value, "Complex<T> requires floating point types.");

    private:
        T real_;
        T imag_;

    public:
        constexpr Complex() noexcept : real_(T(0)), imag_(T(0)) {}
        constexpr Complex(T r) noexcept : real_(r), imag_(T(0)) {}
        constexpr Complex(T r, T i) noexcept : real_(r), imag_(i) {}
        constexpr Complex(const Complex& other) noexcept : real_(other.real_), imag_(other.imag_) {}
        constexpr T real() const noexcept { return real_; }
        constexpr T imag() const noexcept { return imag_; }
        
        void real(T r) noexcept { real_ = r; }
        void imag(T i) noexcept { imag_ = i; }


        Complex operator+(const Complex& other) const noexcept {
            return Complex(real_ + other.real_, imag_ + other.imag_);
        }

        Complex operator-(const Complex& other) const noexcept {
            return Complex(real_ - other.real_, imag_ - other.imag_);
        }

        Complex operator-() const noexcept {
            return Complex(-real_, -imag_);
        }

        Complex operator*(const Complex& other) const noexcept {
            return Complex(
                real_ * other.real_ - imag_ * other.imag_,
                real_ * other.imag_ + imag_ * other.real_
            );
        }

        Complex operator*(T scalar) const noexcept {
            return Complex(real_ * scalar, imag_ * scalar);
        }


        Complex operator/(const Complex& other) const {
            T denom, r, i;
            if (std::abs(other.real_) >= std::abs(other.imag_)) {
                T ratio = other.imag_ / other.real_;
                denom = other.real_ + ratio * other.imag_;
                r = (real_ + imag_ * ratio) / denom;
                i = (imag_ - real_ * ratio) / denom;
            } else {
                T ratio = other.real_ / other.imag_;
                denom = other.imag_ + ratio * other.real_;
                r = (real_ * ratio + imag_) / denom;
                i = (imag_ * ratio - real_) / denom;
            }
            return Complex(r, i);
        }

        Complex operator/(T scalar) const {
            return Complex(real_ / scalar, imag_ / scalar);
        }

        Complex& operator+=(const Complex& other) noexcept {
            real_ += other.real_; imag_ += other.imag_; return *this;
        }
        Complex& operator-=(const Complex& other) noexcept {
            real_ -= other.real_; imag_ -= other.imag_; return *this;
        }
        Complex& operator*=(const Complex& other) noexcept {
            T old_real = real_;
            real_ = old_real * other.real_ - imag_ * other.imag_;
            imag_ = old_real * other.imag_ + imag_ * other.real_;
            return *this;
        }



        Complex conj() const noexcept {
            return Complex(real_, -imag_);
        }

        T normSq() const noexcept {
            return real_ * real_ + imag_ * imag_;
        }

        T abs() const {
            return std::hypot(real_, imag_);
        }

        T arg() const {
            return std::atan2(imag_, real_);
        }

        static Complex exp(const Complex& z) {
            T exp_real = std::exp(z.real_);
            return Complex(exp_real * std::cos(z.imag_), exp_real * std::sin(z.imag_));
        }

        static Complex sqrt(const Complex& z) {
            if (z.real_ == T(0) && z.imag_ == T(0)) return Complex(T(0), T(0));
            T mag = z.abs();
            T r = std::sqrt((mag + z.real_) * T(0.5));
            T i = std::sqrt((mag - z.real_) * T(0.5));
            if (z.imag_ < T(0)) i = -i;
            return Complex(r, i);
        }

        static Complex log(const Complex& z) {
            return Complex(std::log(z.abs()), z.arg());
        }

        static Complex pow(const Complex& base, const Complex& exponent) {
            if (base.abs() < Aurelia::Config::NUMERICAL_EPSILON) return Complex(T(0));
            return exp(exponent * log(base));
        }

        static Complex pow(const Complex& base, T exponent) {
            return pow(base, Complex(exponent));
        }

        bool operator==(const Complex& other) const {
            return std::abs(real_ - other.real_) < Aurelia::Config::NUMERICAL_EPSILON &&
                   std::abs(imag_ - other.imag_) < Aurelia::Config::NUMERICAL_EPSILON;
        }

        friend std::ostream& operator<<(std::ostream& os, const Complex& z) {
            os << "(" << z.real_ << (z.imag_ >= 0 ? "+" : "") << z.imag_ << "i)";
            return os;
        }
    };

    using ComplexH = Complex<long double>;

    template<typename T>
    Complex<T> operator*(T scalar, const Complex<T>& z) {
        return Complex<T>(scalar * z.real(), scalar * z.imag());
    }

} 
} 

#endif