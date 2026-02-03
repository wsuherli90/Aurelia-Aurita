#ifndef AURELIA_ACTIVEMATTER_THERMODYNAMICS_NON_EQUILIBRIUM_H
#define AURELIA_ACTIVEMATTER_THERMODYNAMICS_NON_EQUILIBRIUM_H

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "../../../Foundation/LinearAlgebra/Matrix.h"
#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace ActiveMatter {
namespace Thermodynamics {

    using Real = long double;
    using Matrix3 = Aurelia::Math::Matrix<Real>;

    class DissipationFunction {
    public:
        static Real computeLocalDissipation(const Matrix3& metric_rate, 
                                            const Matrix3& g_inv,
                                            Real viscosity = Aurelia::Config::MESOGLEA_VISCOSITY) {
            
            if (metric_rate.rows() != g_inv.rows())
                throw std::invalid_argument("Dissipation: Metric dimension mismatch.");
            Real geometric_dissipation = 0.0L;
            size_t n = metric_rate.rows(); 
            for(size_t i=0; i<n; ++i) {
                for(size_t j=0; j<n; ++j) {
                    Real rate_mixed_ij = 0.0L;
                    for(size_t k=0; k<n; ++k) {
                        rate_mixed_ij += g_inv(i, k) * metric_rate(k, j);
                    }
                
                    for(size_t k=0; k<n; ++k) {
                        for(size_t l=0; l<n; ++l) {
                            Real dg_dt_ij = metric_rate(i, j);
                            Real dg_dt_kl = metric_rate(k, l);
                            geometric_dissipation += g_inv(i, k) * g_inv(j, l) * dg_dt_ij * dg_dt_kl;
                        }
                    }
                }
            }


            return viscosity * std::max(0.0L, geometric_dissipation);
        }


        static bool checkSecondLaw(Real total_dissipation) {
            return total_dissipation >= -Aurelia::Config::NUMERICAL_EPSILON;
        }

        static Real computeEfficiency(Real work_done, Real chemical_input) {
            if (std::abs(chemical_input) < Aurelia::Config::NUMERICAL_EPSILON) return 0.0L;
            return work_done / chemical_input;
        }
    };

} 
} 
} 

#endif 