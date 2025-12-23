

#ifndef AURELIA_STATMECH_MICROSTRUCTURE_ODF_H
#define AURELIA_STATMECH_MICROSTRUCTURE_ODF_H

#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm> 

#include "../../../Foundation/LinearAlgebra/Vector.h"
#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace StatMech {
namespace Microstructure {

    using Real = long double;
    using Vector3 = Aurelia::Math::Vector<Real>;

 
    struct FiberPopulation {
        Vector3 mu;   
        Real kappa;   
        Real weight;  
        Real normC;   

        FiberPopulation(Vector3 direction, Real k, Real w) 
            : mu(direction.normalized()), kappa(k), weight(w) {
            

            if (kappa < 1.0e-4L) {
                normC = 1.0L / (4.0L * Aurelia::Config::PI);
            } else {
                normC = kappa / (4.0L * Aurelia::Config::PI * std::sinh(kappa));
            }
        }
    };

    class OrientationDistributionFunction {
    private:
        std::vector<FiberPopulation> populations_;

    public:
        OrientationDistributionFunction() {
            populations_.reserve(4); 
        }

 
        void addPopulation(const Vector3& direction, Real alignment_strength, Real fraction) {
            populations_.emplace_back(direction, alignment_strength, fraction);
        }


        void clear() {
            populations_.clear();
        }


        Real evaluate(const Vector3& n) const {
            if (populations_.empty()) {

                return 1.0L / (4.0L * Aurelia::Config::PI);
            }

            Real prob = 0.0L;
            Real total_weight = 0.0L;

            for (const auto& pop : populations_) {

                Real dot_prod = pop.mu.dot(n);
                

                Real pdf = pop.normC * std::exp(pop.kappa * dot_prod);
                
                prob += pop.weight * pdf;
                total_weight += pop.weight;
            }


            if (total_weight > Aurelia::Config::NUMERICAL_EPSILON) {
                return prob / total_weight;
            }
            return prob; 
        }

   
        Real evaluateSpherical(Real theta, Real phi) const {
            Real sin_t = std::sin(theta);
            Vector3 n = {
                sin_t * std::cos(phi),
                sin_t * std::sin(phi),
                std::cos(theta)
            };
            return evaluate(n);
        }

        size_t numPopulations() const { return populations_.size(); }
    };

} 
} 
} 

#endif 