#ifndef LEBESGUE_CONSTANT_H
#define LEBESGUE_CONSTANT_H

#include <cmath>

namespace Aurelia {
namespace Crypto {
namespace Approximation {
namespace Theory {


class LebesgueConstant {
public:
    static double estimate(int n) {

        
        return std::pow(2.0, n / 10.0);
    }
};

} // namespace Theory
} // namespace Approximation
} // namespace Crypto
} // namespace Aurelia

#endif
