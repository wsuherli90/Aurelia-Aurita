

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

            auto integrand = [&](Real theta, Real phi) -> Real {

                Real sin_t = std::sin(theta);
                Vector3 n = {
                    sin_t * std::cos(phi),
                    sin_t * std::sin(phi),
                    std::cos(theta)
                };

                Real projection = std::abs(y.dot(n));
                Real extension = projection * microscopic_scale_;


                Real energy = wlc_model_.computeFreeEnergy(extension);

                return odf_structure_.evaluate(n) * energy;
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