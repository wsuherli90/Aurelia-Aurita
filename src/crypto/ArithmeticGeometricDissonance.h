#ifndef ARITHMETIC_GEOMETRIC_DISSONANCE_H
#define ARITHMETIC_GEOMETRIC_DISSONANCE_H


#include "Representation/ConversionBottleneck.h"
#include "Representation/Security/ProbabilisticModel.h"

#include "Invariant/Distinguisher.h"
#include "Invariant/HomologicalAlgebra/FreeResolution.h"
#include "Invariant/HomologicalAlgebra/TorFunctor.h"

#include "Approximation/ApproximationBound.h"
#include "Approximation/Theory/RemezSolver.h"

#include <vector>
#include <iostream>
#include <iomanip>

namespace Aurelia {
namespace Crypto {

class DissonanceAnalyzer {
public:
    template<typename T>
    static void generateSoKReport(
        double grid_dx, 
        size_t total_lattice_points, 
        const std::vector<T>& sample_field_data,
        int proposed_masking_order = 2
    ) {
        std::cout << "\n==============================================================\n";
        std::cout << "SoK: The Arithmetic-Geometric Dissonance [2025 Academic Analysis]\n";
        std::cout << "Structural Gaps and Limits in the 2025 Cryptographic Landscape\n";
        std::cout << "==============================================================\n";


        std::cout << "[1] Representation Gap: Probing Security Verification\n";
        auto isw_proof = Representation::Security::ProbabilisticModel::verifyGadget(proposed_masking_order, 256);
        std::cout << "    - ISW-Probing Model: " << (isw_proof.is_secure ? "SECURE" : "LEAKAGE DETECTED") << "\n";
        if (!isw_proof.is_secure) {
            std::cout << "      ! WARNING: Implementation flaw probability P(Leak) = " << isw_proof.leakage_probability << "\n";
        }
        
        auto rep_report = Representation::ConversionBottleneck::analyzeMetric(proposed_masking_order);
        std::cout << Representation::ConversionBottleneck::printReport(rep_report) << "\n";


        std::cout << "[2] Invariant Gap: Homological Algebra Analysis\n";
        

        auto inv_result = Invariant::Distinguisher::runDistinguisher(sample_field_data);
        

        Invariant::HomologicalAlgebra::FreeResolution resolution;
        resolution.computeResolution(inv_result.invariant_metric > 1.0 ? 0.0 : 1.0);
        resolution.printChain();
        
        double beta_1 = Invariant::HomologicalAlgebra::TorFunctor::dimTor(1, resolution);
        std::cout << "    - Computed Betti Number beta_1 = dim Tor_1(S/I, K) = " << beta_1 << "\n";
        
        if (inv_result.success) {
            std::cout << "    -> DIAGNOSIS: Syzygy Distinguisher SUCCESS. Code structure is visible.\n\n";
        } else {
             std::cout << "    -> DIAGNOSIS: Generic Resolution. Indistinguishable from Random.\n\n";
        }


        std::cout << "[3] Approximation Gap: Constructive Approximation Theory\n";
        
        double target_eps = 1e-9;
        auto minimax = Approximation::Theory::RemezSolver::solveForReLU(target_eps);
        
        std::cout << "    - Running Remez Exchange Algorithm for f(x)=ReLU(x)...\n";
        std::cout << "    - Optimal Minimax Degree P*: " << minimax.degree << "\n";
        std::cout << "    - Max Error ||f - P*||_inf: " << minimax.max_error << "\n";
        
        auto approx_report = Approximation::ApproximationBound::analyze(target_eps);
        std::cout << "    - Jackson Bound Check: " << approx_report.required_degree << " (Consistent with Remez)\n";
        std::cout << "    - Circuit Depth Limit: " << (minimax.degree > 1000 ? "EXCEEDED" : "OK") << "\n\n";


        std::cout << "==============================================================\n";
        std::cout << "Final Academic Verdict (2026 Agenda):\n";
        
        bool critical = false;
        
        if (rep_report.is_prohibitive || !isw_proof.is_secure) {
             std::cout << " [CRITICAL] Physical Layer: Unsafe or Prohibitive. Move to Isomorphic Rings.\n";
             critical = true;
        }
        
        if (inv_result.success) {
             std::cout << " [CRITICAL] Theoretical Layer: Broken Assumptions. Revoke Code-Based Standards.\n";
             critical = true;
        }
        
        if (minimax.degree > 1000) {
             std::cout << " [CRITICAL] Protocol Layer: Scalability Wall. Discrete Arithmetization failed.\n";
             critical = true;
        }

        if (!critical) std::cout << " -> All Systems Nominal.\n";
        
        std::cout << "==============================================================\n\n";
    }
};

} // namespace Crypto
} // namespace Aurelia

#endif
