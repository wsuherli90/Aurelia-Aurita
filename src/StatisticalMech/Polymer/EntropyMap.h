#ifndef AURELIA_STATMECH_ENTROPY_MAP_H
#define AURELIA_STATMECH_ENTROPY_MAP_H

#include <cmath>
#include <vector>
#include <mutex>

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
    
    struct CachedSphericalNode {
        Vector3 n;   
        Real weight; 
    };

    class EntropyMap {
    private:
        WormLikeChain wlc_model_; 
        Aurelia::StatMech::Microstructure::OrientationDistributionFunction odf_structure_;  
        Real microscopic_scale_;

        static std::vector<CachedSphericalNode> integration_cache_;
        static std::once_flag cache_init_flag_;

    public:
        EntropyMap() : wlc_model_(), odf_structure_() {
            microscopic_scale_ = 0.1L * Aurelia::Config::COLLAGEN_CONTOUR_LENGTH;
            std::call_once(cache_init_flag_, &EntropyMap::initializeCache);
        }

        static void initializeCache() {
            auto nodes = Aurelia::Math::Calculus::GaussLegendre::getSphericalNodes();
            
            integration_cache_.reserve(nodes.size());
            for (const auto& node : nodes) {
                Real theta = node.theta;
                Real phi = node.phi;
                Real w = node.weight;

                Real sin_t = std::sin(theta);
                Vector3 dir = {
                    sin_t * std::cos(phi),
                    sin_t * std::sin(phi),
                    std::cos(theta)
                };
                
                integration_cache_.push_back({dir, w});
            }
        }

        Real operator()(const PointTM& u) const {
            return computeFinslerMetric(u);
        }

        Real computeFinslerMetric(const PointTM& u) const {
            Vector3 y = u.y();
            Real y_mag = y.norm();
            if (y_mag < Aurelia::Config::NUMERICAL_EPSILON) return y_mag;

            Real total_energy = 0.0L;
            size_t cache_size = integration_cache_.size();

            #pragma omp simd reduction(+:total_energy)
            for (size_t i = 0; i < cache_size; ++i) {
                const auto& node = integration_cache_[i];
                const Vector3& n = node.n;

                Real projection = std::abs(y.dot(n));
                Real extension = projection * microscopic_scale_;

                Real wlc_E = wlc_model_.computeFreeEnergy(extension);
                
                Real prob = odf_structure_.evaluate(n);

                total_energy += prob * wlc_E * node.weight;
            }

            return std::sqrt(std::max(0.0L, 2.0L * total_energy));
        }

        void addFiberPopulation(Vector3 axis, Real kappa, Real weight) {
            odf_structure_.addPopulation(axis, kappa, weight);
        }
        void clearStructure() { odf_structure_.clear(); }
    };

    std::vector<CachedSphericalNode> EntropyMap::integration_cache_;
    std::once_flag EntropyMap::cache_init_flag_;

} 
} 
}

#endif