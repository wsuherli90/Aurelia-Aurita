#ifndef AURELIA_ACTIVEMATTER_MORPHOGENESIS_SYMMETRIZATION_H
#define AURELIA_ACTIVEMATTER_MORPHOGENESIS_SYMMETRIZATION_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>

#include "RicciFlow.h"
#include "../../Geometry/Manifold/SlitTangentBundle.h"
#include "../../Foundation/LinearAlgebra/Vector.h"
#include "../../Foundation/Calculus/NumericalDiff.h" 
#include "../../config/BiophyiscalConstants.h"

namespace Aurelia {
namespace ActiveMatter {
namespace Morphogenesis {

    using Real = long double;
    using Vector3 = Aurelia::Math::Vector<Real>;
    using Matrix3 = Aurelia::Math::Matrix<Real>;
    using DiffEngine = Aurelia::Math::Calculus::NumericalDiff;

    class Symmetrization {
    private:
        RicciFlow flow_engine_;
        ChemPot chem_pot_;
        Real symmetry_tolerance_;
        Real tissue_mobility_;

    public:

        Symmetrization() 
            : flow_engine_(), chem_pot_(), 
              symmetry_tolerance_(Aurelia::Config::SYMMETRY_TOLERANCE),
              tissue_mobility_(Aurelia::Config::TISSUE_MOBILITY) {}


        Real computeSymmetryDefect(const std::vector<Real>& curvature_samples) {
            if (curvature_samples.empty()) return 0.0L;

            Real mean = 0.0L;
            for (Real val : curvature_samples) mean += val;
            mean /= curvature_samples.size();

            Real variance = 0.0L;
            for (Real val : curvature_samples) {
                variance += (val - mean) * (val - mean);
            }
            return variance / curvature_samples.size();
        }


        Vector3 computeGradientMu(PointTM& u, Real R_scalar, const Vector3& pos_probe) {
            auto potential_func = [&](const Vector3& p) -> Real {
                PointTM temp = u; temp.set_x(p);
                return chem_pot_.compute(temp, R_scalar);
            };

            Vector3 grad(3);
            for(size_t i=0; i<3; ++i) {
                grad[i] = DiffEngine::diff_1st(potential_func, pos_probe, i);
            }
            return grad;
        }

 
        Vector3 computeVelocity(PointTM& u, Real R_scalar) {
            Vector3 grad_mu = computeGradientMu(u, R_scalar, u.x());
            return grad_mu * (-tissue_mobility_);
        }


        Matrix3 computeVelocityGradient(PointTM& u, Real R_scalar) {
            Matrix3 Dv(3, 3);
            Vector3 x_orig = u.x();
            auto velocity_comp_func = [&](const Vector3& pos, size_t component_idx) -> Real {
                Vector3 grad_at_pos = computeGradientMu(u, R_scalar, pos);
                return grad_at_pos[component_idx] * (-tissue_mobility_);
            };

            for(size_t i=0; i<3; ++i) { 
                auto specific_v_comp = [&](const Vector3& p) { return velocity_comp_func(p, i); };
                
                for(size_t j=0; j<3; ++j) { 
                    Dv(i, j) = DiffEngine::diff_1st(specific_v_comp, x_orig, j);
                }
            }
            return Dv;
        }

        void executeRepair(std::vector<PointTM>& mesh_points,
                           Aurelia::Geometry::Manifold::MetricTensor& g_engine,
                           ChernConn& conn) {
            
            size_t max_steps = Aurelia::Config::MAX_SOLVER_ITERATIONS;
            Real dt = Aurelia::Config::TIME_DT;

            std::cout << "Starting Active Symmetrization..." << std::endl;

            for (size_t t = 0; t < max_steps; ++t) {
                std::vector<Real> curvature_profile;
                std::vector<Vector3> displacements;

                for (size_t i = 0; i < mesh_points.size(); ++i) {
                    PointTM& u = mesh_points[i];

                    g_engine.compute(u);
                    conn.compute(g_engine, u);
                    Matrix3 R_ij = RicciFlow::computeRicciTensor(u, conn, g_engine);
                    Matrix3 g_inv = g_engine.contravariant();
                    Real R_scalar = (g_inv * R_ij).trace(); 
                    
                    curvature_profile.push_back(R_scalar);
                    Vector3 v_flow = computeVelocity(u, R_scalar);
                    Matrix3 Dv_flow = computeVelocityGradient(u, R_scalar);
                    Matrix3 g_new = flow_engine_.evolve(u, g_engine, conn, chem_pot_, v_flow, Dv_flow);
                    displacements.push_back(v_flow * dt);
                }

                for (size_t i = 0; i < mesh_points.size(); ++i) {
                    PointTM& u = mesh_points[i];
                    Vector3 new_pos = u.x() + displacements[i];
                    u.set_x(new_pos);
                }

                Real defect_score = computeSymmetryDefect(curvature_profile);
                
                if (t % 10 == 0) {
                    std::cout << "Time " << t 
                              << " | Sym Defect: " << defect_score 
                              << " | Max v: " << displacements[0].norm()/dt 
                              << std::endl;
                }

                if (defect_score < symmetry_tolerance_) {
                    std::cout << "CONVERGENCE: Soliton State Reached." << std::endl;
                    break;
                }
            }
        }
    };

} 
} 
} 

#endif