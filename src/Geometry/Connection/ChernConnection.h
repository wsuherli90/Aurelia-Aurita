#ifndef AURELIA_GEOMETRY_CHERN_CONNECTION_H
#define AURELIA_GEOMETRY_CHERN_CONNECTION_H

#include <vector>
#include <functional>
#include <cmath>

#include "../Manifold/MetricTensor.h"
#include "../Manifold/SlitTangentBundle.h"
#include "../../Foundation/Calculus/NumericalDiff.h"
#include "../../Foundation/LinearAlgebra/Matrix.h"
#include "../../Foundation/LinearAlgebra/Vector.h"

namespace Aurelia {
namespace Geometry {
namespace Connection {

    using Real = long double;
    using Vector3 = Aurelia::Math::Vector<Real, 3>;
    using Matrix3 = Aurelia::Math::Matrix<Real, 3, 3>;
    using Tensor3 = std::vector<Matrix3>; 
    using Diff = Aurelia::Math::Calculus::NumericalDiff;

    class ChernConnection {
    private:
        Tensor3 gamma_; 
        Matrix3 non_linear_N_;
        Vector3 spray_G_;

        bool computed_;
        Tensor3 computeOsculatingGamma(Aurelia::Geometry::Manifold::MetricTensor& g_engine, 
                                       Aurelia::Geometry::Manifold::PointTM& u) {
            
            const Real h = 1.0e-4L;
            const Real inv_2h = 1.0L / (2.0L * h);

            g_engine.compute(u);
            Matrix3 g_inv = g_engine.contravariant();
            Vector3 x_orig = u.x();


            Tensor3 dg_dx(3, Matrix3(3, 3)); 

            for (size_t k = 0; k < 3; ++k) {

                Vector3 x_p = x_orig; x_p[k] += h; u.set_x(x_p);
                g_engine.compute(u); Matrix3 g_p = g_engine.covariant();

                Vector3 x_m = x_orig; x_m[k] -= h; u.set_x(x_m);
                g_engine.compute(u); Matrix3 g_m = g_engine.covariant();

                dg_dx[k] = (g_p - g_m) * inv_2h;
            }
            u.set_x(x_orig); 
            g_engine.compute(u); 
            Tensor3 osc_gamma(3, Matrix3(3, 3));
            for(size_t i=0; i<3; ++i) {
                for(size_t j=0; j<3; ++j) {
                    for(size_t k=0; k<3; ++k) {
                        Real sum = 0.0L;
                        for(size_t m=0; m<3; ++m) {
                            sum += g_inv(i, m) * (dg_dx[k](m, j) + dg_dx[j](m, k) - dg_dx[m](j, k));
                        }
                        osc_gamma[i](j, k) = 0.5L * sum;
                    }
                }
            }
            return osc_gamma;
        }

    public:
        ChernConnection() : gamma_(3, Matrix3(3, 3)), computed_(false) {}
        void compute(Aurelia::Geometry::Manifold::MetricTensor& g_engine, 
                     Aurelia::Geometry::Manifold::PointTM& u) {
            
            const Real h = 1.0e-4L;
            const Real inv_2h = 1.0L / (2.0L * h);
            Vector3 y_orig = u.y();
            Tensor3 center_gamma = computeOsculatingGamma(g_engine, u);
            
            spray_G_ = Vector3(3, 0.0L);
            for(size_t i=0; i<3; ++i) {
                for(size_t j=0; j<3; ++j) {
                    for(size_t k=0; k<3; ++k) {
                        spray_G_[i] += 0.5L * center_gamma[i](j, k) * y_orig[j] * y_orig[k];
                    }
                }
            }
            auto getSprayVector = [&](Vector3 y_probe) -> Vector3 {

                Aurelia::Geometry::Manifold::PointTM temp_u = u;
                temp_u.set_y(y_probe); 

                Tensor3 temp_gamma = computeOsculatingGamma(g_engine, temp_u);
                
                Vector3 G(3, 0.0L);
                for(size_t i=0; i<3; ++i) {
                    for(size_t j=0; j<3; ++j) {
                        for(size_t k=0; k<3; ++k) {
                            G[i] += 0.5L * temp_gamma[i](j, k) * y_probe[j] * y_probe[k];
                        }
                    }
                }
                return G;
            };

            non_linear_N_ = Matrix3(3, 3, 0.0L);
            for(size_t j=0; j<3; ++j) { 
                Vector3 y_p = y_orig; y_p[j] += h;
                Vector3 G_p = getSprayVector(y_p);

                Vector3 y_m = y_orig; y_m[j] -= h;
                Vector3 G_m = getSprayVector(y_m);

                Vector3 dG_dyj = (G_p - G_m) * inv_2h;

                for(size_t i=0; i<3; ++i) {
                    non_linear_N_(i, j) = dG_dyj[i];
                }
            }
            

            g_engine.compute(u);
            Tensor3 dg_dy(3, Matrix3(3, 3));
            for(size_t m=0; m<3; ++m) {
                Vector3 y_p = y_orig; y_p[m] += h; u.set_y(y_p);
                g_engine.compute(u); Matrix3 g_p = g_engine.covariant();

                Vector3 y_m = y_orig; y_m[m] -= h; u.set_y(y_m);
                g_engine.compute(u); Matrix3 g_m = g_engine.covariant();

                dg_dy[m] = (g_p - g_m) * inv_2h;
            }
            u.set_y(y_orig); 
            g_engine.compute(u); 

            Matrix3 g_inv = g_engine.contravariant();
            Vector3 x_orig = u.x();
            Tensor3 dg_dx(3, Matrix3(3, 3)); 
            for (size_t k = 0; k < 3; ++k) {
                Vector3 x_p = x_orig; x_p[k] += h; u.set_x(x_p);
                g_engine.compute(u); Matrix3 g_p = g_engine.covariant();
                Vector3 x_m = x_orig; x_m[k] -= h; u.set_x(x_m);
                g_engine.compute(u); Matrix3 g_m = g_engine.covariant();
                dg_dx[k] = (g_p - g_m) * inv_2h;
            }
            u.set_x(x_orig);
            g_engine.compute(u);


            Tensor3 delta_g(3, Matrix3(3, 3));
            for(size_t k=0; k<3; ++k) {
                for(size_t i=0; i<3; ++i) {
                    for(size_t j=0; j<3; ++j) {
                        Real d_dx = dg_dx[k](i, j);
                        Real N_correction = 0.0L;
                        for(size_t m=0; m<3; ++m) {
                            N_correction += non_linear_N_(m, k) * dg_dy[m](i, j);
                        }
                        delta_g[k](i, j) = d_dx - N_correction;
                    }
                }
            }


            for(size_t i=0; i<3; ++i) {
                for(size_t j=0; j<3; ++j) {
                    for(size_t k=0; k<3; ++k) {
                        Real sum = 0.0L;
                        for(size_t m=0; m<3; ++m) {
                            sum += g_inv(i, m) * (delta_g[k](m, j) + delta_g[j](m, k) - delta_g[m](j, k));
                        }
                        gamma_[i](j, k) = 0.5L * sum;
                    }
                }
            }

            computed_ = true;
        }


        Real operator()(size_t i, size_t j, size_t k) const {
            if(!computed_) throw std::runtime_error("ChernConnection: Compute not called.");
            return gamma_[i](j, k);
        }


        Real getN(size_t i, size_t j) const {
            if(!computed_) throw std::runtime_error("ChernConnection: Compute not called.");
            return non_linear_N_(i, j);
        }
    };

}
}
} 

#endif