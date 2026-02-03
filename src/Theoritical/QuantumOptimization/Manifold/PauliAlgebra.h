#ifndef AURELIA_THEORETICAL_QUANTUM_PAULI_ALGEBRA_H
#define AURELIA_THEORETICAL_QUANTUM_PAULI_ALGEBRA_H

#include <vector>
#include <string>
#include <complex>
#include <iostream>
#include <array>

namespace Aurelia {
namespace Theoretical {
namespace QuantumOptimization {
namespace Manifold {

    using Complex = std::complex<long double>;

    enum PauliBase { I=0, X=1, Y=2, Z=3 };

    /**
     * @brief PauliString
     * Represents a tensor product of Pauli operators: P = \bigotimes P_i
     * Used to define the Hamiltonian H = \sum c_i O_i and Unitary Generators P_k.
     */
    struct PauliString {
        std::vector<PauliBase> param_string;
        Complex phase; // Global phase (usually 1, i, -1, -i)

        PauliString(size_t num_qubits) : param_string(num_qubits, I), phase(1.0, 0.0) {}

        // Commutator [A, B] = AB - BA
        // For Pauli strings, [P_a, P_b] is either 0 or 2 * P_c
        static std::pair<Complex, PauliString> Commutator(const PauliString& A, const PauliString& B) {
            if (A.param_string.size() != B.param_string.size()) 
                throw std::runtime_error("Pauli dimensions mismatch");

            PauliString res(A.param_string.size());
            Complex phase_accum = A.phase * B.phase;
            
            // Interaction logic
            // X*Y=iZ, Y*X=-iZ, etc.
            // Commutator of individual paulis determines global outcome.
            // If they commute an even number of times, total commutes. Odd -> anticommutes.
            
            int anticommute_count = 0;
            
            for(size_t i=0; i<A.param_string.size(); ++i) {
                PauliBase a = A.param_string[i];
                PauliBase b = B.param_string[i];
                
                if (a == I) res.param_string[i] = b;
                else if (b == I) res.param_string[i] = a;
                else if (a == b) res.param_string[i] = I;
                else {
                    // Non-trivial product
                    anticommute_count++;
                    if (a == X && b == Y) { res.param_string[i] = Z; phase_accum *= Complex(0, 1); }
                    else if (a == Y && b == X) { res.param_string[i] = Z; phase_accum *= Complex(0, -1); }
                    else if (a == Y && b == Z) { res.param_string[i] = X; phase_accum *= Complex(0, 1); }
                    else if (a == Z && b == Y) { res.param_string[i] = X; phase_accum *= Complex(0, -1); }
                    else if (a == Z && b == X) { res.param_string[i] = Y; phase_accum *= Complex(0, 1); }
                    else if (a == X && b == Z) { res.param_string[i] = Y; phase_accum *= Complex(0, -1); }
                }
            }
            
            // [A, B] = AB - BA
            // If A,B commute (anticommute_count even): AB = BA -> [A,B] = 0
            // If A,B anticommute (anticommute_count odd): AB = -BA -> [A,B] = 2AB
            
            if (anticommute_count % 2 == 0) {
                return { Complex(0,0), res };
            } else {
                return { 2.0L * phase_accum, res };
            }
        }
    };

}
}
}
}

#endif // AURELIA_THEORETICAL_QUANTUM_PAULI_ALGEBRA_H
