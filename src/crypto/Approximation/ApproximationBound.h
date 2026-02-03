#ifndef APPROXIMATION_BOUND_H
#define APPROXIMATION_BOUND_H

#include "Manifold.h"
#include "PolynomialSolver.h"
#include <sstream>

namespace Aurelia {
namespace Crypto {
namespace Approximation {


class ApproximationBound {
public:
    struct BoundReport {
        std::string manifold_type;
        double target_epsilon;
        int required_degree;
        double circuit_depth;
        double polynomial_overhead;
    };

    static BoundReport analyze(double epsilon) {
        BoundReport report;
        report.manifold_type = Manifold::getName();
        report.target_epsilon = epsilon;
        
        int k = Manifold::getSmoothnessClass();
        report.required_degree = PolynomialSolver::estimateRequiredDegree(epsilon, k);
        report.circuit_depth = PolynomialSolver::calculateDepth(report.required_degree);
        

        double ideal_d = std::log(1.0/epsilon);
        report.polynomial_overhead = (double)report.required_degree / ideal_d;
        
        return report;
    }

    static std::string report(const BoundReport& r) {
        std::stringstream ss;
        ss << "   [Approximation Bound] Manifold: " << r.manifold_type << "\n";
        ss << "    - Target Epsilon: " << r.target_epsilon << "\n";
        ss << "    - Required Degree: " << r.required_degree << " (Constraint)\n";
        ss << "    - Est. Circ Depth: " << r.circuit_depth << " (Bootstraps required)\n";
        ss << "    - Complexity Penalty: " << r.polynomial_overhead << "x vs Smooth Functions";
        return ss.str();
    }
};

} // namespace Approximation
} // namespace Crypto
} // namespace Aurelia

#endif
