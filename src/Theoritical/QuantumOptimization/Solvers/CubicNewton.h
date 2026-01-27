#ifndef AURELIA_THEORETICAL_QUANTUM_CUBIC_NEWTON_H
#define AURELIA_THEORETICAL_QUANTUM_CUBIC_NEWTON_H

#include "../Sampling/StaticEnvelope.h"
#include "LanczosSolver.h"
#include "../../Manifold/CircuitGraph.h"

namespace Aurelia {
namespace Theoretical {
namespace QuantumOptimization {
namespace Solvers {

    using namespace Aurelia::Math;

    class CubicNewton {
    public:
        using GradientCaller = std::function<std::pair<Real, DynamicVector>(const DynamicVector&)>;
        using JacobianSampler = std::function<DynamicVector(const DynamicVector&, const DynamicVector&, const std::vector<size_t>&, const std::vector<Real>&)>;

    private:
        Sampling::StaticEnvelope envelope_;
        LanczosSolver lanczos_;
        Real L_lipchitz_;
        size_t sample_size_;

    public:
        CubicNewton(Real L = 1.0, size_t sample_s = 200) 
            : L_lipchitz_(L), sample_size_(sample_s), lanczos_(20) {}

        void initialize_envelope(const Manifold::CircuitGraph& graph, const std::vector<std::vector<int>>& measure_qubits) {
            std::cout << "[CubicNewton] Pre-computing Static Envelope...\n";
            envelope_.build_distribution(graph, measure_qubits);
        }

        DynamicVector minimize(DynamicVector x_init, GradientCaller gradOp, JacobianSampler jacobOp) {
            DynamicVector x = x_init;
            
            for(int k=0; k<100; ++k) {
                 auto [fx, gx] = gradOp(x);
                 Real g_norm = gx.norm();

                 std::cout << "Iter K=" << k << " f(x)=" << fx << " |g|=" << g_norm << "\n";
                 
                 if (g_norm < 1e-5) {
                     std::cout << "Converged.\n";
                     break;
                 }

                 std::vector<size_t> indices = envelope_.oblivious_sample(sample_size_);
                 
                 std::vector<Real> weights(indices.size());
                 for(size_t i=0; i<indices.size(); ++i) {
                     Real p = envelope_.get_probability(indices[i]);
                     weights[i] = (p > 1e-12) ? (1.0 / (sample_size_ * p)) : 0.0;
                 }

                 HessianOracle Hv = [&](const DynamicVector& v) {
                     return jacobOp(x, v, indices, weights);
                 };

                 DynamicVector d_star = lanczos_.solve(gx, Hv, L_lipchitz_);

                 x = x + d_star;
            }
            return x;
        }
    };

}
}
}
}

#endif
