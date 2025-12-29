

#ifndef AURELIA_STATMECH_ENTROPY_MAP_H
#define AURELIA_STATMECH_ENTROPY_MAP_H

#include <cmath>
#include <functional>
#include <algorithm>

#include "WormLikeChain.h"
#include "../Microstructure/ODF.h"
#include "../../Foundation/Calculus/Integration.h"
#include "../../Foundation/LinearAlgebra/Vector.h"
#include "../../Geometry/Manifold/SlitTangentBundle.h"
#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace StatMech {
namespace Polymer {

    using Real = long double;
    using Vector3 = Aurelia::Math::Vector<Real>;
    using PointTM = Aurelia::Geometry::Manifold::PointTM;

    using ODFClass = Aurelia::StatMech::Microstructure::OrientationDistributionFunction;

    class EntropyMap {
    private:
        WormLikeChain wlc_model_; 
        ODFClass odf_structure_;  

        Real microscopic_scale_;

    public:
    
        EntropyMap() : wlc_model_(), odf_structure_() {

            microscopic_scale_ = 0.1L * Aurelia::Config::COLLAGEN_CONTOUR_LENGTH;
        }

  
        Real computeFinslerMetric(const PointTM& u) {

            Vector3 y = u.y();
            Real y_mag = y.norm();

            if (y_mag < 1.0e-9L) return y_mag;

            // Cache y components to avoid repeated vector access in integrand
            const Real y0 = y[0], y1 = y[1], y2 = y[2];
            const Real scale = microscopic_scale_;

            auto integrand = [&](Real theta, Real phi) -> Real {
                // Compute direction components without allocating Vector3
                Real sin_t = std::sin(theta);
                Real cos_t = std::cos(theta);
                Real cos_p = std::cos(phi);
                Real sin_p = std::sin(phi);
                
                Real nx = sin_t * cos_p;
                Real ny = sin_t * sin_p;
                Real nz = cos_t;

                // Inline dot product computation
                Real projection = std::abs(y0 * nx + y1 * ny + y2 * nz);
                Real extension = projection * scale;

                Real energy = wlc_model_.computeFreeEnergy(extension);

                // Use optimized ODF evaluate without Vector allocation
                return odf_structure_.evaluate(nx, ny, nz) * energy;
            };

            Real total_energy = Aurelia::Math::Calculus::GaussLegendre::integrateSpherical(integrand);


            return std::sqrt(std::max(0.0L, 2.0L * total_energy));
        }
        

        void addFiberPopulation(Vector3 axis, Real kappa, Real weight) {
            odf_structure_.addPopulation(axis, kappa, weight);
        }

 
        void clearStructure() {
            odf_structure_.clear();
        }

    
        void setMicroscopicScale(Real scale_ratio) {
            microscopic_scale_ = scale_ratio * Aurelia::Config::COLLAGEN_CONTOUR_LENGTH;
        }
        

        const WormLikeChain& getPolymerModel() const {
            return wlc_model_;
        }
    };

} 
} 
}

#endif 