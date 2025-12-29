

#ifndef AURELIA_STATMECH_WORMLIKECHAIN_H
#define AURELIA_STATMECH_WORMLIKECHAIN_H

#include <cmath>
#include <stdexcept>
#include <algorithm> 

#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace StatMech {
namespace Polymer {

    using Real = long double;

    class WormLikeChain {
    private:
        Real Lp_; 
        Real Lc_; 
        Real T_; 
        Real kBT_; 

    public:
     
        WormLikeChain(Real Lp = Aurelia::Config::COLLAGEN_PERSISTENCE_LENGTH,
                      Real Lc = Aurelia::Config::COLLAGEN_CONTOUR_LENGTH,
                      Real T  = Aurelia::Config::TEMPERATURE_K)
            : Lp_(Lp), Lc_(Lc), T_(T) {
            
            kBT_ = Aurelia::Config::BOLTZMANN_CONSTANT * T_;
        }

    
        Real computeForce(Real z) const noexcept {

            if (z >= Lc_ * Aurelia::Config::WLC_SINGULARITY_CUTOFF) 
                return Aurelia::Config::WLC_FORCE_PENALTY; 
            
            if (z <= 0.0L) return 0.0L; // Slack chain

            Real ratio = z / Lc_;
            Real term1 = 1.0L / (4.0L * (1.0L - ratio) * (1.0L - ratio));
            Real term2 = -0.25L;
            Real term3 = ratio;

            return (kBT_ / Lp_) * (term1 + term2 + term3);
        }

        Real computeFreeEnergy(Real z) const noexcept {

            if (z >= Lc_ * Aurelia::Config::WLC_SINGULARITY_CUTOFF) 
                return Aurelia::Config::WLC_ENERGY_BARRIER; 
            
            if (z <= 0.0L) return 0.0L;

            Real ratio = z / Lc_; 
            
  
            Real term1 = 1.0L / (4.0L * (1.0L - ratio));
            

            Real term2 = -ratio / 4.0L;
            

            Real term3 = (ratio * ratio) / 2.0L;

            Real constant = -0.25L;

            Real prefactor = (kBT_ * Lc_) / Lp_;

            return prefactor * (term1 + term2 + term3 + constant);
        }

        // Getters for properties
        Real stiffness() const { return kBT_ / Lp_; }
        Real limit() const { return Lc_; }
        Real thermalEnergy() const { return kBT_; }
    };

} 
} 
} 

#endif 