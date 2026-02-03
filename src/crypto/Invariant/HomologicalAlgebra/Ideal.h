#ifndef IDEAL_H
#define IDEAL_H

#include "PolynomialRing.h"
#include <vector>

namespace Aurelia {
namespace Crypto {
namespace Invariant {
namespace HomologicalAlgebra {


class Ideal {
public:
    std::vector<Polynomial> generators;
    bool is_groebner_basis = false;

    void addGenerator(const Polynomial& p) {
        generators.push_back(p);
        is_groebner_basis = false;
    }

    void computeGroebnerBasis() {

        

        
        is_groebner_basis = true;
    }
    
    int getBasisSize() const {
        return generators.size();
    }
};

} // namespace HomologicalAlgebra
} // namespace Invariant
} // namespace Crypto
} // namespace Aurelia

#endif
