#ifndef MANIFOLD_H
#define MANIFOLD_H

#include <cmath>
#include <string>

namespace Aurelia {
namespace Crypto {
namespace Approximation {


class Manifold {
public:
    static double evaluateReLU(double x) {
        return (x > 0.0) ? x : 0.0;
    }


    static int getSmoothnessClass() {
        return 1;
    }

    static std::string getName() {
        return "ReLU Manifold (Non-Smooth)";
    }
};

} // namespace Approximation
} // namespace Crypto
} // namespace Aurelia

#endif
