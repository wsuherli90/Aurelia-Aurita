#ifndef AURELIA_FOUNDATION_DYNAMIC_VECTOR_H
#define AURELIA_FOUNDATION_DYNAMIC_VECTOR_H

#include <vector>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <initializer_list>
#include <algorithm>

namespace Aurelia {
namespace Math {

    using Real = long double;

    class DynamicVector {
    private:
        std::vector<Real> data_;

    public:
        DynamicVector() = default;
        explicit DynamicVector(size_t size, Real val = 0.0) : data_(size, val) {}
        DynamicVector(std::vector<Real> d) : data_(std::move(d)) {}
        DynamicVector(std::initializer_list<Real> list) : data_(list) {}

        static DynamicVector Zero(size_t n) {
            return DynamicVector(n, 0.0);
        }

        Real& operator[](size_t i) { return data_[i]; }
        const Real& operator[](size_t i) const { return data_[i]; }
        size_t size() const { return data_.size(); }

        DynamicVector operator+(const DynamicVector& other) const {
            if(size() != other.size()) throw std::invalid_argument("Vector size mismatch (+)");
            DynamicVector res(size());
            for(size_t i=0; i<size(); ++i) res[i] = data_[i] + other[i];
            return res;
        }

        DynamicVector operator-(const DynamicVector& other) const {
            if(size() != other.size()) throw std::invalid_argument("Vector size mismatch (-)");
            DynamicVector res(size());
            for(size_t i=0; i<size(); ++i) res[i] = data_[i] - other[i];
            return res;
        }

        DynamicVector operator-() const {
            DynamicVector res(size());
            for(size_t i=0; i<size(); ++i) res[i] = -data_[i];
            return res;
        }

        DynamicVector operator*(Real scalar) const {
            DynamicVector res(size());
            for(size_t i=0; i<size(); ++i) res[i] = data_[i] * scalar;
            return res;
        }

        friend DynamicVector operator*(Real scalar, const DynamicVector& v) {
            return v * scalar;
        }

        Real dot(const DynamicVector& other) const {
            if(size() != other.size()) throw std::invalid_argument("Vector size mismatch (dot)");
            Real sum = 0.0;
            for(size_t i=0; i<size(); ++i) sum += data_[i] * other[i];
            return sum;
        }

        Real squaredNorm() const {
            Real sum = 0.0;
            for(auto x : data_) sum += x * x;
            return sum;
        }

        Real norm() const {
            return std::sqrt(squaredNorm());
        }
        
        friend std::ostream& operator<<(std::ostream& os, const DynamicVector& v) {
            os << "[";
            for(size_t i=0; i<v.size(); ++i) {
                os << v[i] << (i < v.size()-1 ? ", " : "");
            }
            os << "]";
            return os;
        }
    };

}
}

#endif
