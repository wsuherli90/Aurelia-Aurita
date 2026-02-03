#ifndef POLYNOMIAL_RING_H
#define POLYNOMIAL_RING_H

#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <iostream>

namespace Aurelia {
namespace Crypto {
namespace Invariant {
namespace HomologicalAlgebra {


struct Monomial {
    std::vector<int> exponents;
    
    bool operator<(const Monomial& other) const {
        int deg_a = 0, deg_b = 0;
        for(int e : exponents) deg_a += e;
        for(int e : other.exponents) deg_b += e;
        
        if (deg_a != deg_b) return deg_a < deg_b;
        

        for (int i = exponents.size() - 1; i >= 0; --i) {
            if (exponents[i] != other.exponents[i]) return exponents[i] > other.exponents[i];
        }
        return false;
    }
};

class Polynomial {
public:
    std::map<Monomial, double> terms;

    void addTerm(const Monomial& m, double coeff) {
        if (std::abs(coeff) < 1e-9) return;
        terms[m] += coeff;
        if (std::abs(terms[m]) < 1e-9) terms.erase(m);
    }
    
    std::string toString() const {
        if (terms.empty()) return "0";
        std::string s = "";
        for (auto it = terms.rbegin(); it != terms.rend(); ++it) {
            s += std::to_string((int)it->second) + "*X^" +  std::to_string(it->first.exponents[0]) + " ";
        }
        return s;
    }
    
    bool isZero() const { return terms.empty(); }
};

} // namespace HomologicalAlgebra
} // namespace Invariant
} // namespace Crypto
} // namespace Aurelia

#endif
