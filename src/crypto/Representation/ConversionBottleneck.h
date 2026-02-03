#ifndef CONVERSION_BOTTLENECK_H
#define CONVERSION_BOTTLENECK_H

#include "BooleanGadget.h"
#include "ArithmeticRing.h"
#include <sstream>
#include <iomanip>

namespace Aurelia {
namespace Crypto {
namespace Representation {


class ConversionBottleneck {
public:
    struct BottleneckReport {
        int masking_order;
        double b2a_cycles;
        double a2b_cycles;
        double total_overhead_factor;
        bool is_prohibitive;
    };

    static BottleneckReport analyzeMetric(int t) {
        BottleneckReport report;
        report.masking_order = t;


        
        double bit_width = ArithmeticRing::getBitWidth();
        double adder_cost = 5.0 * BooleanGadget::estimateAndGateCost(t) + 
                            10.0 * BooleanGadget::estimateXorGateCost(t); 
        
        report.b2a_cycles = bit_width * adder_cost;
        

        report.a2b_cycles = report.b2a_cycles * 1.5;

        double baseline_op = ArithmeticRing::estimateMultiplicationCost();
        
        report.total_overhead_factor = (report.b2a_cycles + report.a2b_cycles) / baseline_op;
        

        report.is_prohibitive = (report.total_overhead_factor > 100.0);

        return report;
    }

    static std::string printReport(const BottleneckReport& r) {
        std::stringstream ss;
        ss << "   [G_B2A Analysis] Masking Order t=" << r.masking_order << "\n";
        ss << "    - Logic Gates Cost: " << std::fixed << std::setprecision(1) << r.b2a_cycles << " cycles (B2A)\n";
        ss << "    - Arithmetic Cost:  " << r.a2b_cycles << " cycles (A2B)\n";
        ss << "    - Overhead Factor:  " << r.total_overhead_factor << "x vs Unmasked\n";
        if (r.is_prohibitive) {
            ss << "    -> STATUS: PROHIBITIVE. Representation Gap is critical barrier.\n";
        } else {
            ss << "    -> STATUS: MANAGEABLE. Hardware acceleration feasible.\n";
        }
        return ss.str();
    }
};

} // namespace Representation
} // namespace Crypto
} // namespace Aurelia

#endif
