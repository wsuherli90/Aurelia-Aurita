#ifndef AURELIA_STATMECH_WLC_H
#define AURELIA_STATMECH_WLC_H

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace StatMech {
namespace Polymer {

    using Real = long double;

    class WormLikeChain {
    private:
        Real persistence_length_;
        Real contour_length_;
        Real temperature_;
        Real kbT_;

    public:
        
        WormLikeChain() {
            persistence_length_ = Aurelia::Config::COLLAGEN_PERSISTENCE_LENGTH;
            contour_length_     = Aurelia::Config::COLLAGEN_CONTOUR_LENGTH;
            temperature_        = Aurelia::Config::PHYSIOLOGICAL_TEMP;
            kbT_                = Aurelia::Config::BOLTZMANN_CONST * temperature_;
        }

        #pragma omp declare simd uniform(this) linear(extension)
        Real computeFreeEnergy(Real extension) const {
            
            Real x = std::abs(extension);
            
            Real L = contour_length_;
            if (x >= L * 0.99L) x = L * 0.99L; 

            Real ratio = x / L;
            Real term1 = 1.0L - ratio;
            Real term1_sq = term1 * term1;

            Real factor = kbT_ * contour_length_ / (4.0L * persistence_length_);
            
            Real partA = 1.0L / term1_sq;
            
            
            Real energy = (kbT_ / persistence_length_) * (0.25L * partA - 0.25L + ratio);

            return energy;
        }

        void setPersistenceLength(Real Lp) { persistence_length_ = Lp; }
    };

} 
} 
}

#endif