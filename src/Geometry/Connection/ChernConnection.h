#ifndef AURELIA_GEOMETRY_CHERN_CONNECTION_H
#define AURELIA_GEOMETRY_CHERN_CONNECTION_H

#include <vector>
#include <cmath>
#include <array>

#include "../Manifold/MetricTensor.h"
#include "../../Foundation/LinearAlgebra/Matrix.h"
#include "../../Foundation/LinearAlgebra/Vector.h"
#include "../../Foundation/LinearAlgebra/Tensor4.h"
#include "../../Foundation/Calculus/StencilCache.h"

namespace Aurelia {
namespace Geometry {
namespace Connection {

    using Real = long double;
    using Vector3 = Aurelia::Math::Vector<Real, 3>;
    using Matrix3 = Aurelia::Math::Matrix<Real, 3, 3>;
    using Tensor3 = std::vector<Matrix3>; 
    using Stencil = Aurelia::Math::Calculus::MetricStencilCache;


    class ChernConnection {
    private:
        Tensor3 gamma_;         
        Matrix3 non_linear_N_;  
        bool computed_;

        Tensor3 computeOsculatingGamma(const Stencil& stencil, const Matrix3& g_inv) {
            Tensor3 osc_gamma(3, Matrix3(3, 3));
            
            for(size_t k=0; k<3; ++k) {
                for(size_t i=0; i<3; ++i) {
                    for(size_t j=0; j<3; ++j) {
    
                        Real sum = 0.0L;
                        for(size_t m=0; m<3; ++m) {
                            Real d_j_gim = stencil.partialDerivative(i, m, j);
                            Real d_i_gjm = stencil.partialDerivative(j, m, i);
                            Real d_m_gij = stencil.partialDerivative(i, j, m);
                            
                            sum += g_inv(k, m) * (d_j_gim + d_i_gjm - d_m_gij);
                        }
                        osc_gamma[k](i, j) = 0.5L * sum;
                    }
                }
            }
            return osc_gamma;
        }


        Vector3 computeSpray(const Tensor3& osc_gamma, const Vector3& y) {
            Vector3 G(0.0L);
            for(size_t i=0; i<3; ++i) {
                for(size_t j=0; j<3; ++j) {
                    for(size_t k=0; k<3; ++k) {
                        G[i] += 0.5L * osc_gamma[i](j, k) * y[j] * y[k];
                    }
                }
            }
            return G;
        }

    public:
        ChernConnection() : gamma_(3, Matrix3(3, 3)), computed_(false) {}

        template <typename MetricEngine>
        void compute(MetricEngine& g_engine, Aurelia::Geometry::Manifold::PointTM& u) {
            
            Real h = Aurelia::Config::FINITE_DIFFERENCE_STEP;
            Vector3 y_orig = u.y();
            
            Stencil center_stencil;
            center_stencil.populate(g_engine, u, Stencil::SPATIAL_X);
            Matrix3 g_center = center_stencil.getCenter();
            Matrix3 g_inv = g_center.inverse();
            
            non_linear_N_ = Matrix3(0.0L);

            const Real inv_12h = 1.0L / (12.0L * h);

            Tensor3 dg_dy(3, Matrix3(3,3)); 

            auto computePerturbed = [&](const Vector3& y_p) -> std::pair<Vector3, Matrix3> {
                PointTM u_p = u; u_p.set_y(y_p);
                Stencil stencil_p;
                stencil_p.populate(g_engine, u_p, Stencil::SPATIAL_X);
                Matrix3 g_p = stencil_p.getCenter();
                Matrix3 g_inv_p = g_p.inverse();
                Tensor3 gamma_p = computeOsculatingGamma(stencil_p, g_inv_p);
                Vector3 G_p = computeSpray(gamma_p, y_p);
                return {G_p, g_p};
            };

            for(size_t j=0; j<3; ++j) { 
                
                Vector3 y_p2 = y_orig; y_p2[j] += 2.0L*h;
                auto [G_p2, g_p2] = computePerturbed(y_p2);

                Vector3 y_p1 = y_orig; y_p1[j] += h;
                auto [G_p1, g_p1] = computePerturbed(y_p1);

                Vector3 y_m1 = y_orig; y_m1[j] -= h;
                auto [G_m1, g_m1] = computePerturbed(y_m1);

                Vector3 y_m2 = y_orig; y_m2[j] -= 2.0L*h;
                auto [G_m2, g_m2] = computePerturbed(y_m2);

                Vector3 dG_dy = (-G_p2 + G_p1 * 8.0L - G_m1 * 8.0L + G_m2) * inv_12h;
                
                for(size_t i=0; i<3; ++i) {
                    non_linear_N_(i, j) = dG_dy[i];
                }

                for(size_t r=0; r<3; ++r) {
                    for(size_t c=0; c<3; ++c) {
                        Real val = (-g_p2(r,c) + 8.0L*g_p1(r,c) - 8.0L*g_m1(r,c) + g_m2(r,c)) * inv_12h;
                        dg_dy[j](r, c) = val;
                    }
                }
            }
            
            u.set_y(y_orig);

            Tensor3 delta_g_dx(3, Matrix3(3,3)); 

            for(size_t k=0; k<3; ++k) { 
                for(size_t i=0; i<3; ++i) {
                    for(size_t j=0; j<3; ++j) {
                        
                        Real partial_x = center_stencil.partialDerivative(i, j, k);
                        
                        Real correction = 0.0L;
                        for(size_t m=0; m<3; ++m) {
                            Real partial_y = dg_dy[m](i, j);
                            correction += non_linear_N_(m, k) * partial_y;
                        }

                        delta_g_dx[k](i, j) = partial_x - correction;
                    }
                }
            }

            for(size_t k=0; k<3; ++k) {
                for(size_t i=0; i<3; ++i) {
                    for(size_t j=0; j<3; ++j) {
                        Real sum = 0.0L;
                        for(size_t m=0; m<3; ++m) {
                            Real term1 = delta_g_dx[i](j, m);
                            Real term2 = delta_g_dx[j](i, m);
                            Real term3 = delta_g_dx[m](i, j);
                            
                            sum += g_inv(k, m) * (term1 + term2 - term3);
                        }
                        gamma_[k](i, j) = 0.5L * sum;
                    }
                }
            }

            computed_ = true;
        }

        Real operator()(size_t k, size_t i, size_t j) const {
            #ifdef DEBUG_MATH
            if(!computed_) throw std::runtime_error("ChernConnection: Compute not called.");
            #endif
            return gamma_[k](i, j);
        }

        Real getN(size_t i, size_t j) const {
             #ifdef DEBUG_MATH
            if(!computed_) throw std::runtime_error("ChernConnection: Compute not called.");
            #endif
            return non_linear_N_(i, j);
        }
    };

}
}
} 

#endif