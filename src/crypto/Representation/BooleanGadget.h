#ifndef BOOLEAN_GADGET_H
#define BOOLEAN_GADGET_H

#include <cmath>
#include <vector>
#include <string>

namespace Aurelia {
namespace Crypto {
namespace Representation {


class BooleanGadget {
public:
    static double estimateAndGateCost(int masking_order_t) {

        double base_cycles = 1.0;
        double quadratic_overhead = std::pow(masking_order_t + 1, 2.0);
        return base_cycles * quadratic_overhead * 2.5;
    }

    static double estimateXorGateCost(int masking_order_t) {

        return (double)(masking_order_t + 1);
    }
    
    static std::string getGadgetType() { return "ISW-2003 (Proven Secure)"; }
};

} // namespace Representation
} // namespace Crypto
} // namespace Aurelia

#endif
