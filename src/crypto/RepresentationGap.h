#ifndef REPRESENTATION_GAP_H
#define REPRESENTATION_GAP_H

#include <cmath>
#include <string>
#include <sstream>
#include <iostream>

namespace Aurelia {
namespace Crypto {


class RepresentationGap {
public:
    struct SecurityMetrics {
        int masking_order_t;
        double conversion_overhead_factor;
        unsigned long long estimated_cycles;
        std::string security_level;
        std::string barrier_description;
    };

    static constexpr double ML_KEM_Q = 3329.0;
    static constexpr double LOG_Q = 11.70099;
    static constexpr double CHES_GADGET_CONSTANT = 12.5; 

    static SecurityMetrics calculateGap(int masking_order_t, double modulus_q_override, size_t operation_count) {
        SecurityMetrics metrics;
        metrics.masking_order_t = masking_order_t;


        double t_factor = static_cast<double>(masking_order_t * masking_order_t);
        double structural_penalty = t_factor * LOG_Q;
        

        metrics.conversion_overhead_factor = 1.0 + (structural_penalty * CHES_GADGET_CONSTANT);
        
        metrics.estimated_cycles = static_cast<unsigned long long>(operation_count * metrics.conversion_overhead_factor);

        if (metrics.conversion_overhead_factor < 15.0) {
            metrics.security_level = "Low Assurance (Leaky)";
            metrics.barrier_description = "Manageable software overhead.";
        } else if (metrics.conversion_overhead_factor < 150.0) {
            metrics.security_level = "Standard Assurance (Masked)";
            metrics.barrier_description = "Significant arithmetic bottleneck.";
        } else {
            metrics.security_level = "High Assurance (Prohibitive)";
            metrics.barrier_description = "Representation Gap exceeds hardware limits.";
        }

        return metrics;
    }

    static std::string report(const SecurityMetrics& m) {
        std::stringstream ss;
        ss << "[Representation Gap] Physical Layer (Boolean-Arithmetic Dissonance):\n";
        ss << "   Masking Order (t):      " << m.masking_order_t << "\n";
        ss << "   Cost Model:             O(t^2 * log q)\n";
        ss << "   Overhead Factor:        " << m.conversion_overhead_factor << "x\n";
        ss << "   Est. Cycles:            " << m.estimated_cycles << "\n";
        ss << "   Status:                 " << m.security_level << "\n";
        ss << "   Conclusion:             " << m.barrier_description;
        return ss.str();
    }
};

} 
} 

#endif 
