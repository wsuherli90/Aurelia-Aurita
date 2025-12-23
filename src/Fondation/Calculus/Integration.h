#ifndef AURELIA_FOUNDATION_INTEGRATION_H
#define AURELIA_FOUNDATION_INTEGRATION_H

#include <vector>
#include <functional>
#include <cmath>

#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace Math {
namespace Calculus {

    using Real = long double;

    class GaussLegendre {
    private:
        static constexpr size_t N_POINTS = 10;
        
        static constexpr Real NODES[N_POINTS] = {
            -0.9739065285171717L, 0.9739065285171717L,
            -0.8650633666889845L, 0.8650633666889845L,
            -0.6794095682990244L, 0.6794095682990244L,
            -0.4333953941292472L, 0.4333953941292472L,
            -0.1488743389816312L, 0.1488743389816312L
        };

        static constexpr Real WEIGHTS[N_POINTS] = {
            0.0666713443086881L, 0.0666713443086881L,
            0.1494513491505806L, 0.1494513491505806L,
            0.2190863625159820L, 0.2190863625159820L,
            0.2692667193099963L, 0.2692667193099963L,
            0.2955242247147529L, 0.2955242247147529L
        };

    public:
        static Real integrate(const std::function<Real(Real)>& func, Real a, Real b) {
            Real half_len = (b - a) * 0.5L;
            Real center = (b + a) * 0.5L;
            Real sum = 0.0L;

            for (size_t i = 0; i < N_POINTS; ++i) {
                Real x_real = half_len * NODES[i] + center;
                sum += WEIGHTS[i] * func(x_real);
            }

            return half_len * sum;
        }
        static Real integrateSpherical(const std::function<Real(Real, Real)>& func) {
            Real theta_a = 0.0L, theta_b = Aurelia::Config::PI;
            Real phi_a = 0.0L, phi_b = Aurelia::Config::TWO_PI;

            Real theta_half = (theta_b - theta_a) * 0.5L;
            Real theta_center = (theta_b + theta_a) * 0.5L;
            
            Real phi_half = (phi_b - phi_a) * 0.5L;
            Real phi_center = (phi_b + phi_a) * 0.5L;

            Real sum = 0.0L;

            for (size_t i = 0; i < N_POINTS; ++i) {
                Real theta = theta_half * NODES[i] + theta_center;
                Real w_theta = WEIGHTS[i];
                Real jacobian = std::sin(theta); 

                for (size_t j = 0; j < N_POINTS; ++j) {
                    Real phi = phi_half * NODES[j] + phi_center;
                    Real w_phi = WEIGHTS[j];

                    sum += w_theta * w_phi * func(theta, phi) * jacobian;
                }
            }

            return theta_half * phi_half * sum;
        }
    };

} 
} 
} 

#endif 