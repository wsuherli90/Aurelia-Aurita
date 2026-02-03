#ifndef TOR_FUNCTOR_H
#define TOR_FUNCTOR_H

#include "FreeResolution.h"

namespace Aurelia {
namespace Crypto {
namespace Invariant {
namespace HomologicalAlgebra {


class TorFunctor {
public:
    static double dimTor(int i, const FreeResolution& res) {
        if (i < 0 || i >= res.chain.size()) return 0.0;
        
        return (double)res.chain[i].rank;
    }
    
    static double computeEulerCharacteristic(const FreeResolution& res) {
        double chi = 0.0;
        for (size_t i = 0; i < res.chain.size(); ++i) {
            double sign = (i % 2 == 0) ? 1.0 : -1.0;
            chi += sign * res.chain[i].rank;
        }
        return chi;
    }
};

} // namespace HomologicalAlgebra
} // namespace Invariant
} // namespace Crypto
} // namespace Aurelia

#endif
