#ifndef AURELIA_THEORETICAL_QUANTUM_CIRCUIT_GRAPH_H
#define AURELIA_THEORETICAL_QUANTUM_CIRCUIT_GRAPH_H

#include <vector>
#include <set>
#include <iostream>
#include <queue>

namespace Aurelia {
namespace Theoretical {
namespace QuantumOptimization {
namespace Manifold {

    class CircuitGraph {
    public:
        struct GateNode {
            int id;
            std::vector<int> qubits_involved;
            int layer_depth;
        };

    private:
        int num_qubits_;
        std::vector<GateNode> gates_;
        std::vector<std::vector<int>> qubit_to_gates_; 

    public:
        CircuitGraph(int n_qubits) : num_qubits_(n_qubits), qubit_to_gates_(n_qubits) {}

        void addGate(int gate_id, const std::vector<int>& qubits, int layer) {
            gates_.push_back({gate_id, qubits, layer});
            for(int q : qubits) {
                qubit_to_gates_[q].push_back(gate_id);
            }
        }

        int getLightConeSize(const std::vector<int>& measured_qubits) const {
            if (measured_qubits.empty()) return 0;

            std::set<int> visited_gates;
            std::queue<int> q;

            std::vector<bool> qubit_active(num_qubits_, false);
            for(int qb : measured_qubits) {
                if(qb >= 0 && qb < num_qubits_) qubit_active[qb] = true;
            }

            int count = 0;

            for(auto it = gates_.rbegin(); it != gates_.rend(); ++it) {
                const auto& gate = *it;
                
                bool touches_active_lightcone = false;
                for(int q_gate : gate.qubits_involved) {
                    if (qubit_active[q_gate]) {
                        touches_active_lightcone = true;
                        break;
                    }
                }
                
                if (touches_active_lightcone) {
                    count++;
                    
                    for(int q_gate : gate.qubits_involved) {
                        qubit_active[q_gate] = true;
                    }
                }
            }
            return count;
        }
    };

}
}
}
}

#endif
