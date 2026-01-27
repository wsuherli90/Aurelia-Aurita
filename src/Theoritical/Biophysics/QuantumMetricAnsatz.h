#ifndef AURELIA_THEORETICAL_BIOPHYSICS_QUANTUM_METRIC_ANSATZ_H
#define AURELIA_THEORETICAL_BIOPHYSICS_QUANTUM_METRIC_ANSATZ_H

#include <vector>
#include <cmath>
#include <complex>

#include "../../Foundation/LinearAlgebra/DynamicVector.h"
#include "../../ActiveMatter/Morphogenesis/RicciFlow.h"
#include "../QuantumOptimization/Manifold/CircuitGraph.h"

namespace Aurelia {
namespace Theoretical {
namespace Biophysics {

    using Real = long double;
    using Complex = std::complex<Real>;
    using namespace Aurelia::ActiveMatter::Morphogenesis;
    using namespace Aurelia::Theoretical::QuantumOptimization::Manifold;

    class QuantumMetricAnsatz {
    private:
        CircuitGraph circuit_graph_;
        size_t num_qubits_;
        size_t n_params_;
        Real scale_epsilon_;

    public:
        QuantumMetricAnsatz(size_t n_qubits, size_t n_layers) 
            : num_qubits_(n_qubits), 
              circuit_graph_(n_qubits), 
              scale_epsilon_(0.01) 
        {
            int p_idx = 0;
            for(size_t l=0; l<n_layers; ++l) {
                for(size_t q=0; q<n_qubits; ++q) {
                    circuit_graph_.addGate(p_idx++, {static_cast<int>(q)}, static_cast<int>(l));
                }
                for(size_t q=0; q<n_qubits-1; ++q) {
                    circuit_graph_.addGate(-1, {static_cast<int>(q), static_cast<int>(q+1)}, static_cast<int>(l));
                }
            }
            n_params_ = p_idx;
        }

        size_t getNumParams() const { return n_params_; }
        const CircuitGraph& getGraph() const { return circuit_graph_; }

        void decode(const Aurelia::Math::DynamicVector& theta, MetricFieldSoA& field) const {
            
            size_t total_points = field.total_size;
            size_t points_per_qubit = (total_points + num_qubits_ - 1) / num_qubits_;

            #pragma omp parallel for
            for(size_t i=0; i<total_points; ++i) {
                size_t q_idx = i / points_per_qubit;
                if(q_idx >= num_qubits_) q_idx = num_qubits_ - 1;

                Real param = 0.0;
                if(q_idx < theta.size()) param = theta[q_idx]; 

                Real perturbation = scale_epsilon_ * std::cos(param);

                Matrix3 g = field.get(i);
                
                g(0,0) += perturbation;
                g(1,1) += perturbation;
                g(2,2) += perturbation;
                
                field.set(i, g);
            }
        }

        void accumulate_gradient_contribution(
            const Aurelia::Math::DynamicVector& theta, 
            const MetricFieldSoA& grad_wrt_metric, 
            Aurelia::Math::DynamicVector& grad_wrt_theta) const 
        {
            size_t total_points = grad_wrt_metric.total_size;
            size_t points_per_qubit = (total_points + num_qubits_ - 1) / num_qubits_;
            
            std::vector<Real> local_grads(num_qubits_, 0.0);

            #pragma omp parallel for reduction(+:local_grads[:num_qubits_])
            for(size_t i=0; i<total_points; ++i) {
                size_t q_idx = i / points_per_qubit;
                if(q_idx >= num_qubits_) q_idx = num_qubits_ - 1;

                Real param = (q_idx < theta.size()) ? theta[q_idx] : 0.0;
                
                Real d_metric_d_theta = -scale_epsilon_ * std::sin(param);
                
                Matrix3 dL_dg = grad_wrt_metric.get(i); 
                
                Real term = (dL_dg(0,0) + dL_dg(1,1) + dL_dg(2,2)) * d_metric_d_theta;
                
                local_grads[q_idx] += term;
            }
            
            for(size_t q=0; q<num_qubits_; ++q) {
                if(q < grad_wrt_theta.size()) {
                    grad_wrt_theta[q] += local_grads[q];
                }
            }
        }
    };

}
}
}

#endif
