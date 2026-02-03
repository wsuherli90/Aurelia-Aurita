#ifndef INVARIANT_GAP_H
#define INVARIANT_GAP_H

#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>

namespace Aurelia {
namespace Crypto {


class InvariantGap {
public:
    struct BettiMetrics {
        double structural_entropy;
        double syzygy_dissonance;
        double betti_number_beta_1;
        bool distinguisher_success;
    };

    template<typename Vector3Type>
    static BettiMetrics analyzeStructure(const std::vector<Vector3Type>& field_data) {
        BettiMetrics metrics;
        
        if (field_data.empty()) {
            metrics.structural_entropy = 1.0; 
            metrics.syzygy_dissonance = 0.0;
            metrics.betti_number_beta_1 = 0.0;
            metrics.distinguisher_success = false;
            return metrics;
        }


        double M[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        
        for (const auto& vec : field_data) {
            double v[3] = { static_cast<double>(vec[0]), static_cast<double>(vec[1]), static_cast<double>(vec[2]) };
            for (int r = 0; r < 3; r++) {
                for (int c = 0; c < 3; c++) {
                    M[r][c] += v[r] * v[c];
                }
            }
        }

        double norm_factor = 1.0 / field_data.size();
        for(int r=0; r<3; r++) for(int c=0; c<3; c++) M[r][c] *= norm_factor;


        double evals[3];
        calculateEigenvaluesDirect(M, evals);
        
        // Sort eigenvalues descending
        std::sort(evals, evals+3, std::greater<double>());
        double lambda1 = evals[0];


        
        double trace = M[0][0] + M[1][1] + M[2][2];
        
        if (std::abs(trace) > 1e-12) {

            double order_param = 1.5 * (lambda1 - (1.0/3.0));
            if (order_param < 0) order_param = 0;
            if (order_param > 1) order_param = 1;

            metrics.syzygy_dissonance = order_param;
            metrics.betti_number_beta_1 = lambda1;
        } else {
            metrics.syzygy_dissonance = 0.0;
            metrics.betti_number_beta_1 = 0.0;
        }

        metrics.structural_entropy = 1.0 - metrics.syzygy_dissonance;
        

        metrics.distinguisher_success = (metrics.syzygy_dissonance > 0.6);

        return metrics;
    }

    static std::string report(const BettiMetrics& m) {
        std::stringstream ss;
        ss << "[Invariant Gap] Theoretical Layer (Syzygy Distinguisher):\n";
        ss << "   Betti Number (beta_1):  " << std::fixed << std::setprecision(4) << m.betti_number_beta_1 << "\n";
        ss << "   Syzygy Dissonance:      " << m.syzygy_dissonance << " (Order Parameter)\n";
        ss << "   Structural Variance:    " << m.structural_entropy << "\n";
        ss << "   Distinguisher Result:   " << (m.distinguisher_success ? "STRUCTURED (Distinguishable)" : "GENERIC (Indistinguishable)");
        return ss.str();
    }

    template<typename Vector3>
    static double calculateLocalEntropy(const std::vector<Vector3>& directors) {
        BettiMetrics m = analyzeStructure(directors);
        return m.structural_entropy; 
    }

private:
    static void calculateEigenvaluesDirect(const double A[3][3], double out_lambda[3]) {
        double p1 = A[0][1]*A[0][1] + A[0][2]*A[0][2] + A[1][2]*A[1][2];
        if (p1 < 1e-16) {

            out_lambda[0] = A[0][0];
            out_lambda[1] = A[1][1];
            out_lambda[2] = A[2][2];
            return;
        }

        double q = (A[0][0] + A[1][1] + A[2][2]) / 3.0;
        double p2 = (A[0][0] - q)*(A[0][0] - q) + (A[1][1] - q)*(A[1][1] - q) + (A[2][2] - q)*(A[2][2] - q) + 2.0 * p1;
        double p = std::sqrt(p2 / 6.0);
        

        double B[3][3];
        double inv_p = 1.0 / p;
        for(int r=0; r<3; r++) {
            for(int c=0; c<3; c++) {
                B[r][c] = (r==c ? (A[r][c] - q) : A[r][c]) * inv_p;
            }
        }

        double detB = B[0][0] * (B[1][1] * B[2][2] - B[1][2] * B[2][1]) -
                      B[0][1] * (B[1][0] * B[2][2] - B[1][2] * B[2][0]) +
                      B[0][2] * (B[1][0] * B[2][1] - B[1][1] * B[2][0]);

        double r = detB / 2.0;


        double phi = 0.0;
        if (r <= -1) phi = 3.14159265358979323846 / 3.0;
        else if (r >= 1) phi = 0.0;
        else phi = std::acos(r) / 3.0;

        out_lambda[0] = q + 2.0 * p * std::cos(phi);
        out_lambda[1] = q + 2.0 * p * std::cos(phi + 2.0 * 3.14159265358979323846 / 3.0);
        out_lambda[2] = q + 2.0 * p * std::cos(phi + 4.0 * 3.14159265358979323846 / 3.0);


    }
};

} 
} 

#endif 
