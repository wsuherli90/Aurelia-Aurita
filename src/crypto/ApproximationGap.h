#ifndef APPROXIMATION_GAP_H
#define APPROXIMATION_GAP_H

#include <cmath>
#include <string>
#include <sstream>
#include <iostream>

namespace Aurelia {
namespace Crypto {


class ApproximationGap {
public:
    struct ApproximationMetrics {
        double grid_resolution_step; 
        double target_epsilon;
        double required_degree_d;
        double actual_precision_error;
        double dissonance_cost;
        std::string scaling_law;
    };

    static ApproximationMetrics calculateDissonance(double dx, int smoothness_k) {
        ApproximationMetrics metrics;
        metrics.grid_resolution_step = dx;
        
        double n = 1.0 / (dx > 0 ? dx : 0.001);
        metrics.actual_precision_error = std::pow(n, -static_cast<double>(smoothness_k));
        

        metrics.target_epsilon = 1e-9;
        

        metrics.required_degree_d = std::pow(metrics.target_epsilon, -1.0 / (double)smoothness_k);
        

        metrics.dissonance_cost = metrics.required_degree_d / n;
        
        if (smoothness_k == 1) metrics.scaling_law = "Linear (Prohibitive)";
        else metrics.scaling_law = "Super-Polynomial (Manageable)";

        return metrics;
    }

    static std::string report(const ApproximationMetrics& m) {
        std::stringstream ss;
        ss << "[Approximation Gap] Protocol Layer (Manifold Mismatch):\n";
        ss << "   Grid Step (dx):         " << m.grid_resolution_step << "\n";
        ss << "   Discrete Error:         " << m.actual_precision_error << " (via Jackson's Thm)\n";
        ss << "   Required Degree (d):    " << m.required_degree_d << " (for upsilon=" << m.target_epsilon << ")\n";
        ss << "   Scaling Law:            " << m.scaling_law << "\n";
        ss << "   Dissonance Cost:        " << m.dissonance_cost << "x (Circuit Depth Penalty)";
        return ss.str();
    }

    static double calculateLocalDissonance(double torsion_magnitude, double dx) {
        double manifold_complexity = 1.0 + torsion_magnitude;
        double scaling_factor = std::pow(manifold_complexity, 0.5); 
        
        return scaling_factor; 
    }
};

} 
} 

#endif 
