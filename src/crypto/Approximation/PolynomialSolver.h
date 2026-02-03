#ifndef POLYNOMIAL_SOLVER_H
#define POLYNOMIAL_SOLVER_H

#include <cmath>

namespace Aurelia {
namespace Crypto {
namespace Approximation {


class PolynomialSolver {
public:
    static constexpr double PI = 3.14159265358979323846;

    static int estimateRequiredDegree(double target_epsilon, int smoothness_k) {

        
        if (target_epsilon <= 0) return 1000000;
        
        double exponent = 1.0 / (double)smoothness_k;
        double d = std::pow(1.0 / target_epsilon, exponent);
        

        return static_cast<int>(std::ceil(d));
    }
    
    static double calculateDepth(int degree) {

        return std::ceil(std::log2((double)degree));
    }
};

} // namespace Approximation
} // namespace Crypto
} // namespace Aurelia

#endif
