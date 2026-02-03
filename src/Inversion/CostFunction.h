#ifndef AURELIA_INVERSION_COST_FUNCTION_H
#define AURELIA_INVERSION_COST_FUNCTION_H

#include <vector>
#include <cmath>
#include <functional>
#include <stdexcept>

#include "../../../Foundation/LinearAlgebra/Vector.h"
#include "../../../Foundation/LinearAlgebra/Matrix.h"
#include "../../../Geometry/Manifold/MetricTensor.h"
#include "../../../ActiveMatter/Thermodynamics/ChemicalPotential.h"
#include "../../../config/BiophyiscalConstants.h"

namespace Aurelia {
namespace Inversion {

    using Real = long double;
    using VectorX = Aurelia::Math::Vector<Real>;
    using Matrix3 = Aurelia::Math::Matrix<Real>;
    using MetricEngine = Aurelia::Geometry::Manifold::MetricTensor;

    class CostFunction {
    private:
        Real lambda_data_;   
        Real lambda_thermo_; 
        Real lambda_smooth_; 

    public:
   
        CostFunction(Real l_data = 1.0L, Real l_thermo = 0.1L, Real l_smooth = 0.05L)
            : lambda_data_(l_data), lambda_thermo_(l_thermo), lambda_smooth_(l_smooth) {}

    
        Real computeDataMisfit(const VectorX& sim_state, 
                               const VectorX& obs_data, 
                               const std::vector<Matrix3>& metric_field) const {
            
            if (sim_state.size() != obs_data.size())
                throw std::invalid_argument("CostFunction: Data dimension mismatch.");

            size_t n_dofs = sim_state.size();
            size_t n_points = n_dofs / 3;

            if (metric_field.size() != n_points)
                throw std::invalid_argument("CostFunction: Metric field resolution mismatch.");

            Real total_sq_error = 0.0L;


            for (size_t k = 0; k < n_points; ++k) {

                size_t idx_x = 3 * k;
                size_t idx_y = 3 * k + 1;
                size_t idx_z = 3 * k + 2;

                Real dx = sim_state[idx_x] - obs_data[idx_x];
                Real dy = sim_state[idx_y] - obs_data[idx_y];
                Real dz = sim_state[idx_z] - obs_data[idx_z];


                const Matrix3& g = metric_field[k];


                Real term_x = g(0,0)*dx + g(0,1)*dy + g(0,2)*dz;
                Real term_y = g(1,0)*dx + g(1,1)*dy + g(1,2)*dz;
                Real term_z = g(2,0)*dx + g(2,1)*dy + g(2,2)*dz;

                Real local_norm_sq = dx*term_x + dy*term_y + dz*term_z;

                total_sq_error += local_norm_sq;
            }

            return 0.5L * total_sq_error;
        }

    
        Real computeThermodynamicPrior(const std::vector<Matrix3>& active_stress_field) const {
            Real penalty = 0.0L;
            

            Real limit = Aurelia::Config::MAX_CELL_STRESS_LIMIT;
            Real limit_sq = limit * limit;

            for (const auto& stress : active_stress_field) {

                Real stress_mag_sq = 0.0L;
                for(size_t i=0; i<3; ++i) {
                    for(size_t j=0; j<3; ++j) {
                        Real s = stress(i,j);
                        stress_mag_sq += s * s;
                    }
                }
                

                penalty += std::exp(stress_mag_sq / limit_sq);
            }
            return penalty;
        }

   
        Real computeGeometricRegularization(const std::vector<Matrix3>& metric_field) const {
            Real roughness = 0.0L;
            Real dx = Aurelia::Config::GRID_DX;
            Real inv_dx_sq = 1.0L / (dx * dx);


            for (size_t i = 1; i < metric_field.size(); ++i) {
                // Difference Matrix
                Matrix3 diff = metric_field[i] - metric_field[i-1];
                

                Real norm_sq = 0.0L;
                for(size_t r=0; r<3; ++r) {
                    for(size_t c=0; c<3; ++c) {
                        Real val = diff(r,c);
                        norm_sq += val * val;
                    }
                }
                

                roughness += norm_sq * inv_dx_sq;
            }

            return roughness;
        }

  
        Real computeTotalCost(const VectorX& sim, 
                              const VectorX& obs, 
                              const std::vector<Matrix3>& active_stress,
                              const std::vector<Matrix3>& metric_field) const {

            Real J1 = computeDataMisfit(sim, obs, metric_field);
            
            Real J2 = computeThermodynamicPrior(active_stress);
            
            Real J3 = computeGeometricRegularization(metric_field);

            return J1 + (lambda_thermo_ * J2) + (lambda_smooth_ * J3);
        }
    };

} 
} 

#endif 