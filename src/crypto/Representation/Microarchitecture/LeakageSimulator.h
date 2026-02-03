#ifndef LEAKAGE_SIMULATOR_H
#define LEAKAGE_SIMULATOR_H

#include <vector>
#include <cmath>

namespace Aurelia {
namespace Crypto {
namespace Representation {
namespace Microarchitecture {


class LeakageSimulator {
public:
    static double calculateSNR(int masking_order) {

        
        double effective_noise = std::pow(masking_order, 2.0);
        return 1.0 / effective_noise;
    }
};

} // namespace Microarchitecture
} // namespace Representation
} // namespace Crypto
} // namespace Aurelia

#endif
