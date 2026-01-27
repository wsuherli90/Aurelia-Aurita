#ifndef AURELIA_THEORETICAL_QUANTUM_STATIC_ENVELOPE_H
#define AURELIA_THEORETICAL_QUANTUM_STATIC_ENVELOPE_H

#include <vector>
#include <random>
#include <numeric>
#include <algorithm>
#include "../Manifold/CircuitGraph.h"

namespace Aurelia {
namespace Theoretical {
namespace QuantumOptimization {
namespace Sampling {

    using Real = long double;

    class StaticEnvelope {
    private:
        std::vector<Real> probabilities_;
        std::vector<Real> alias_probs_; 
        std::vector<int> alias_table_;
        std::mt19937 gen_;
        bool initialized_ = false;

    public:
        StaticEnvelope() : gen_(std::random_device{}()) {}

        Real get_probability(size_t index) const {
            if(index >= probabilities_.size()) return 0.0;
            return probabilities_[index];
        }

        void build_distribution(const Manifold::CircuitGraph& graph, 
                                const std::vector<std::vector<int>>& measurements_qubits) {
            
            size_t M = measurements_qubits.size();
            probabilities_.resize(M);
            
            Real sum = 0.0;
            for(size_t i=0; i<M; ++i) {
                int degree = graph.getLightConeSize(measurements_qubits[i]);
                probabilities_[i] = (degree > 0) ? (Real)degree : 1.0; 
                sum += probabilities_[i];
            }

            for(auto& p : probabilities_) p /= sum;

            prepare_alias_method();
            initialized_ = true;
        }

        std::vector<size_t> oblivious_sample(size_t sample_size) {
            if(!initialized_) throw std::runtime_error("StaticEnvelope not initialized!");
            
            std::vector<size_t> indices(sample_size);
            std::uniform_real_distribution<Real> dist_u(0.0, 1.0);
            std::uniform_int_distribution<size_t> dist_i(0, probabilities_.size() - 1);

            for(size_t i=0; i<sample_size; ++i) {
                size_t k = dist_i(gen_);
                Real u = dist_u(gen_);
                
                if (u < alias_probs_[k]) {
                    indices[i] = k;
                } else {
                    indices[i] = alias_table_[k];
                }
            }
            return indices;
        }

    private:
        void prepare_alias_method() {
            size_t n = probabilities_.size();
            alias_probs_.resize(n);
            alias_table_.resize(n);
            
            std::vector<size_t> small, large;
            for(size_t i=0; i<n; ++i) {
                alias_probs_[i] = probabilities_[i] * n; 
                if(alias_probs_[i] < 1.0) small.push_back(i);
                else large.push_back(i);
            }
            
            while(!small.empty() && !large.empty()) {
                size_t l = small.back(); small.pop_back();
                size_t g = large.back(); large.pop_back();
                
                alias_table_[l] = g;
                alias_probs_[g] = (alias_probs_[g] + alias_probs_[l]) - 1.0;
                
                if(alias_probs_[g] < 1.0) small.push_back(g);
                else large.push_back(g);
            }
             while(!large.empty()) {
                size_t g = large.back(); large.pop_back();
                alias_probs_[g] = 1.0;
            }
            while(!small.empty()) {
                size_t l = small.back(); small.pop_back();
                alias_probs_[l] = 1.0;
            }
        }
    };

}
}
}
}

#endif
