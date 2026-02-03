#ifndef DISTINGUISHER_H
#define DISTINGUISHER_H

#include "SyzygyBasis.h"
#include "BettiTable.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace Aurelia {
namespace Crypto {
namespace Invariant {


class Distinguisher {
public:
    struct AttackResult {
        bool success;
        double confidence;
        double invariant_metric;
        std::string diagnosis;
    };

    template<typename Vector3Type>
    static AttackResult runDistinguisher(const std::vector<Vector3Type>& field_data) {

        double M[3][3] = {{0}};
        for(const auto& v : field_data) {
             double d[3] = { (double)v[0], (double)v[1], (double)v[2] };
             for(int r=0; r<3; r++) for(int c=0; c<3; c++) M[r][c] += d[r]*d[c];
        }
        double norm = 1.0/field_data.size();
        for(int r=0; r<3; r++) for(int c=0; c<3; c++) M[r][c] *= norm;
        
        double evals[3];
        solveEigensystem(M, evals);
        std::sort(evals, evals+3, std::greater<double>());


        auto basis_props = SyzygyBasis::extractFromDirectors({evals[0], evals[1], evals[2]});


        auto betti_stats = BettiTable::computeBettiNumbers(basis_props);


        AttackResult result;
        result.invariant_metric = betti_stats.beta_1_1;
        result.success = (result.invariant_metric > 1.0);
        
        if (result.success) {
            result.confidence = 0.99;
            result.diagnosis = "Anomalous Syzygies Detected! Code is distinguishable from Random.";
        } else {
            result.confidence = 0.05;
            result.diagnosis = "Generic Betti Table. Indistinguishable from Random.";
        }
        
        return result;
    }

private:

    static void solveEigensystem(const double A[3][3], double out_lambda[3]) {
        double p1 = A[0][1]*A[0][1] + A[0][2]*A[0][2] + A[1][2]*A[1][2];
        double q = (A[0][0] + A[1][1] + A[2][2]) / 3.0;
        if (p1 < 1e-16) { out_lambda[0]=A[0][0]; out_lambda[1]=A[1][1]; out_lambda[2]=A[2][2]; return; }
        double p2 = (A[0][0] - q)*(A[0][0] - q) + (A[1][1] - q)*(A[1][1] - q) + (A[2][2] - q)*(A[2][2] - q) + 2.0 * p1;
        double p = std::sqrt(p2 / 6.0);
        double B[3][3]; double inv_p = 1.0/p;
        for(int r=0; r<3; r++) for(int c=0; c<3; c++) B[r][c] = (r==c ? (A[r][c]-q) : A[r][c]) * inv_p;
        double detB = B[0][0] * (B[1][1] * B[2][2] - B[1][2] * B[2][1]) - B[0][1] * (B[1][0] * B[2][2] - B[1][2] * B[2][0]) + B[0][2] * (B[1][0] * B[2][1] - B[1][1] * B[2][0]);
        double phi = std::acos(std::max(-1.0, std::min(1.0, detB/2.0))) / 3.0;
        out_lambda[0] = q + 2.0 * p * std::cos(phi);
        out_lambda[1] = q + 2.0 * p * std::cos(phi + 2.09439510239);
        out_lambda[2] = q + 2.0 * p * std::cos(phi + 4.18879020479);
    }
};

} // namespace Invariant
} // namespace Crypto
} // namespace Aurelia

#endif
