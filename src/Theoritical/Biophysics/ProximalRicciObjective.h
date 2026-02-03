#ifndef AURELIA_THEORETICAL_BIOPHYSICS_PROXIMAL_RICCI_OBJECTIVE_H
#define AURELIA_THEORETICAL_BIOPHYSICS_PROXIMAL_RICCI_OBJECTIVE_H

#include <utility>
#include "../../Foundation/LinearAlgebra/DynamicVector.h"
#include "../../ActiveMatter/Morphogenesis/RicciFlow.h"
#include "../../ActiveMatter/Thermodynamics/ChemicalPotential.h"
#include "QuantumMetricAnsatz.h"

namespace Aurelia {
namespace Theoretical {
namespace Biophysics {

    using Real = long double;
    using DynamicVector = Aurelia::Math::DynamicVector;
    using namespace Aurelia::ActiveMatter::Morphogenesis;
    
    class ProximalRicciObjective {
    private:
        RicciFlow& flow_engine_;
        const QuantumMetricAnsatz& ansatz_;
        MetricFieldSoA& g_current_; 
        const MetricFieldSoA& g_prev_;   
        const Aurelia::ActiveMatter::Thermodynamics::ChemicalPotential& chem_pot_;
        Real dt_;

        MetricFieldSoA rate_buffer_;
        RicciFlow::Workspace workspace_;
        
    public:
        ProximalRicciObjective(
            RicciFlow& engine,
            const QuantumMetricAnsatz& ansatz,
            MetricFieldSoA& current_buffer,
            const MetricFieldSoA& prev_field,
            const Aurelia::ActiveMatter::Thermodynamics::ChemicalPotential& chem,
            Real dt
        ) : flow_engine_(engine), 
            ansatz_(ansatz), 
            g_current_(current_buffer), 
            g_prev_(prev_field), 
            chem_pot_(chem),
            dt_(dt),
            rate_buffer_(current_buffer.nx, current_buffer.ny, current_buffer.nz) 
        {
             workspace_.resize(current_buffer.total_size);
        }

        std::pair<Real, DynamicVector> operator()(const DynamicVector& theta) {
            
            #pragma omp parallel for
            for(size_t i=0; i<g_current_.total_size; ++i) {
                Matrix3 id; id(0,0)=1; id(1,1)=1; id(2,2)=1;
                g_current_.set(i, id);
            }
            ansatz_.decode(theta, g_current_);

            flow_engine_.computeFlowRate(g_current_, rate_buffer_, workspace_, chem_pot_);
            
            MetricFieldSoA grad_g(g_current_.nx, g_current_.ny, g_current_.nz);
            Real loss = 0.0L;
            Real inv_dt = 1.0L / dt_;

            #pragma omp parallel for reduction(+:loss)
            for(size_t i=0; i<g_current_.total_size; ++i) {
                Matrix3 rate = rate_buffer_.get(i); 
                Matrix3 g = g_current_.get(i);
                Matrix3 g0 = g_prev_.get(i);
                
                Matrix3 residual = g - g0 - (rate * dt_);
                
                loss += 0.5L * residual.squaredNorm();
                
                grad_g.set(i, residual);
            }
            
            DynamicVector grad_theta(theta.size(), 0.0);
            ansatz_.accumulate_gradient_contribution(theta, grad_g, grad_theta);
            
            return {loss, grad_theta};
        }
    };
    
    class ProximalRicciJacobian {
    public:
         DynamicVector operator()(const DynamicVector& x, const DynamicVector& v,
                                  const std::vector<size_t>& indices, const std::vector<Real>& weights) {
             return v; 
         }
    };

}
}
}

#endif
