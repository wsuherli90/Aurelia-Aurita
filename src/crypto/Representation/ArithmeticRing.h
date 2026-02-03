#ifndef ARITHMETIC_RING_H
#define ARITHMETIC_RING_H

#include <cmath>

namespace Aurelia {
namespace Crypto {
namespace Representation {


class ArithmeticRing {
public:
    static constexpr double MODULUS_Q = 3329.0;
    static constexpr int N = 256;

    static double estimateMultiplicationCost() {
        return 1.0;
    }

    static double estimateAdditionCost() {
        return 1.0;
    }
    
    static double getBitWidth() {
        return std::log2(MODULUS_Q);
    }
};

} // namespace Representation
} // namespace Crypto
} // namespace Aurelia

#endif
