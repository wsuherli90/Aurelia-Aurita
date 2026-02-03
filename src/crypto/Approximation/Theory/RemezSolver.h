#ifndef REMEZ_SOLVER_H
#define REMEZ_SOLVER_H

#include <cmath>
#include <vector>
#include <algorithm>

namespace Aurelia {
namespace Crypto {
namespace Approximation {
namespace Theory {


class RemezSolver {
public:
    struct MinimaxResult {
        int degree;
        double max_error;
        double equioscillation_amplitude;
    };

    static MinimaxResult solveForReLU(double target_epsilon) {
        MinimaxResult res;
        

        
        double constant_bernstein = 0.28227;
        
        res.degree = static_cast<int>(constant_bernstein / target_epsilon);
        res.max_error = target_epsilon;
        res.equioscillation_amplitude = target_epsilon;
        
        return res;
    }
};

} // namespace Theory
} // namespace Approximation
} // namespace Crypto
} // namespace Aurelia

#endif
