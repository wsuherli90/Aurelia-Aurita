#ifndef AURELIA_ACTIVEMATTER_MORPHOGENESIS_RICCI_FLOW_H
#define AURELIA_ACTIVEMATTER_MORPHOGENESIS_RICCI_FLOW_H

#include <vector>
#include <cmath>

#include "../Thermodynamics/ChemicalPotential.h"
#include "../../Geometry/Manifold/MetricTensor.h"
#include "../../Geometry/Connection/ChernConnection.h"
#include "../../Foundation/Calculus/NumericalDiff.h"
#include "../../Foundation/LinearAlgebra/Vector.h"
#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace ActiveMatter {
namespace Morphogenesis {

    using Real = long double;
    using Vector3 = Aurelia::Math::Vector<Real, 3>;
    using Matrix3 = Aurelia::Math::Matrix<Real, 3, 3>;
    using PointTM = Aurelia::Geometry::Manifold::PointTM;
    using ChernConn = Aurelia::Geometry::Connection::ChernConnection;
    using ChemPot = Aurelia::ActiveMatter::Thermodynamics::ChemicalPotential;

    class RicciFlow {
    private:
        Real dt_; 

    public:
        explicit RicciFlow(Real time_step = Aurelia::Config::TIME_DT) 
            : dt_(time_step) {}
        static Matrix3 computeRicciTensor(PointTM& u, 
                                          ChernConn& conn, 
                                          Aurelia::Geometry::Manifold::MetricTensor& g_engine) {
            
            Matrix3 Ricci(3, 3, 0.0L);
            Real h = Aurelia::Config::FINITE_DIFFERENCE_STEP; 
            auto getGamma = [&](Vector3 x_pos, size_t k_upper, size_t i_lower, size_t j_lower) -> Real {
                PointTM temp = u; temp.set_x(x_pos);
                g_engine.compute(temp);
                conn.compute(g_engine, temp);
                return conn(k_upper, i_lower, j_lower);
            };

            Vector3 x = u.x();


            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = 0; j < 3; ++j) {
                    Real sum = 0.0L;
                    for (size_t k = 0; k < 3; ++k) {
                        Vector3 x_pk = x; x_pk[k] += h;
                        Vector3 x_mk = x; x_mk[k] -= h;
                        Real d_k_Gamma_ij = (getGamma(x_pk, k, i, j) - getGamma(x_mk, k, i, j)) / (2.0L * h);

                        Vector3 x_pj = x; x_pj[j] += h;
                        Vector3 x_mj = x; x_mj[j] -= h;
                        Real d_j_Gamma_ik = (getGamma(x_pj, k, i, k) - getGamma(x_mj, k, i, k)) / (2.0L * h);
                        conn.compute(g_engine, u); 
                        
                        Real term_quad = 0.0L;
                        for(size_t m=0; m<3; ++m) {
                            term_quad += conn(k, k, m) * conn(m, i, j) 
                                       - conn(k, i, m) * conn(m, k, j);
                        }

                        sum += (d_k_Gamma_ij - d_j_Gamma_ik + term_quad);
                    }
                    Ricci(i, j) = sum;
                }
            }
            return Ricci;
        }

        static Matrix3 computeActiveHessian(PointTM& u,
                                            ChernConn& conn,
                                            const ChemPot& chem_pot,
                                            Real R_scalar) {
            Matrix3 H(3, 3, 0.0L);
            Vector3 x = u.x();
            auto getMu = [&](Vector3 pos) -> Real {
                PointTM temp = u; temp.set_x(pos);
                return chem_pot.compute(temp, R_scalar);
            };

            Vector3 d_mu(3);
            for(size_t a=0; a<3; ++a) {
                d_mu[a] = Aurelia::Math::Calculus::NumericalDiff::diff_1st(
                    [&](const Vector3& p){ return getMu(p); }, x, a);
            }

            for(size_t a=0; a<3; ++a) {
                for(size_t b=0; b<3; ++b) {
                    Real d2_mu;
                    if(a == b) 
                        d2_mu = Aurelia::Math::Calculus::NumericalDiff::diff_2nd(
                            [&](const Vector3& p){ return getMu(p); }, x, a);
                    else
                        d2_mu = Aurelia::Math::Calculus::NumericalDiff::diff_mixed(
                            [&](const Vector3& p){ return getMu(p); }, x, a, b);

                    Real conn_correction = 0.0L;
                    for(size_t c=0; c<3; ++c) {
                        conn_correction += conn(c, a, b) * d_mu[c];
                    }

                    H(a, b) = d2_mu - conn_correction;
                }
            }
            return H;
        }

        static Matrix3 computeLieDerivative(const Matrix3& g, 
                                            const Matrix3& D_v) { 
            Matrix3 Lie(3, 3);
            for(size_t a=0; a<3; ++a) {
                for(size_t b=0; b<3; ++b) {
                    Real term1 = 0.0L; 
                    Real term2 = 0.0L; 
                    for(size_t k=0; k<3; ++k) {
                        term1 += g(b, k) * D_v(k, a); 
                        term2 += g(a, k) * D_v(k, b); 
                    }
                    Lie(a, b) = term1 + term2;
                }
            }
            return Lie;
        }


        Matrix3 evolve(PointTM& u, 
                       Aurelia::Geometry::Manifold::MetricTensor& g_engine,
                       ChernConn& conn,
                       const ChemPot& chem_pot,
                       const Vector3& 
                       const Matrix3& D_flow_velocity) { 
           

            // Ensure pointed and optimize 
            Matrix3 g_old = g_engine.covariant();
            Matrix3 g_inv = g_engine.contravariant();
            Matrix3 R_ij = computeRicciTensor(u, conn, g_engine);
            Real R_scalar = (g_inv * R_ij).trace();
            Matrix3 ActiveStress = computeActiveHessian(u, conn, chem_pot, R_scalar);
            Matrix3 LieDrag = computeLieDerivative(g_old, D_flow_velocity);
            Matrix3 rate_of_change = (R_ij * -2.0L) + ActiveStress + LieDrag;
            Matrix3 g_new = g_old + (rate_of_change * dt_);

            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = i + 1; j < 3; ++j) {
                    Real avg = (g_new(i, j) + g_new(j, i)) * 0.5L;
                    g_new(i, j) = avg;
                    g_new(j, i) = avg;
                }
            }

            return g_new;
        }
    };

} 
} 
} 

#endif