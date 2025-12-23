#ifndef AURELIA_GEOMETRY_FLAG_CURVATURE_H
#define AURELIA_GEOMETRY_FLAG_CURVATURE_H

#include <vector>
#include <cmath>
#include <iostream>
#include <functional>

#include "../Connection/ChernConnection.h"
#include "../Manifold/MetricTensor.h"
#include "../../Foundation/Calculus/NumericalDiff.h"
#include "../../Foundation/LinearAlgebra/Vector.h"
#include "../../Foundation/LinearAlgebra/Matrix.h"
#include "../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace Geometry {
namespace Invariants {

    using Real = long double;
    using Vector3 = Aurelia::Math::Vector<Real>;
    using Matrix3 = Aurelia::Math::Matrix<Real>;
    using Diff = Aurelia::Math::Calculus::NumericalDiff;
    
    using namespace Aurelia::Geometry::Connection;
    using namespace Aurelia::Geometry::Manifold;

    class FlagCurvature {
    public:

        static Real compute(MetricTensor& metric, 
                            ChernConnection& conn, 
                            PointTM& u, 
                            const Vector3& V) {
            

            metric.compute(u);
            const Matrix3& g = metric.covariant();

            auto innerProduct = [&](const Vector3& a, const Vector3& b) -> Real {
                Real sum = 0.0L;
                for(size_t i=0; i<3; ++i)
                    for(size_t j=0; j<3; ++j)
                        sum += g(i, j) * a[i] * b[j];
                return sum;
            };

            Real g_yy = innerProduct(u.y(), u.y());
            Real g_VV = innerProduct(V, V);
            Real g_yV = innerProduct(u.y(), V);

            Real denominator = g_yy * g_VV - g_yV * g_yV;

            if (std::abs(denominator) < Aurelia::Config::NUMERICAL_EPSILON) {
                return 0.0L; 
            }


            auto getGammaFunc_X = [&](size_t i, size_t j, size_t l) {
                return [&](const Vector3& x_pos) -> Real {
                    PointTM temp = u; temp.set_x(x_pos);
                    metric.compute(temp);
                    conn.compute(metric, temp);
                    return conn(i, j, l);
                };
            };

            auto getGammaFunc_Y = [&](size_t i, size_t j, size_t l) {
                return [&](const Vector3& y_dir) -> Real {
                    PointTM temp = u; temp.set_y(y_dir);
                    metric.compute(temp);
                    conn.compute(metric, temp);
                    return conn(i, j, l);
                };
            };

            metric.compute(u);
            conn.compute(metric, u);


            Vector3 W(3, 0.0L);
            Vector3 y = u.y();

            
            for (size_t i = 0; i < 3; ++i) { 
                Real sum_R = 0.0L;

                for (size_t j = 0; j < 3; ++j) {
                    for (size_t k = 0; k < 3; ++k) {
                        for (size_t l = 0; l < 3; ++l) {
 
                            Real factor = y[j] * V[k] * y[l];
                            
                            if (std::abs(factor) < 1.0e-15L) continue;

                            Real term_comm = 0.0L;
                            for (size_t m = 0; m < 3; ++m) {
                                term_comm += conn(i, m, k) * conn(m, j, l) 
                                           - conn(i, m, l) * conn(m, j, k);
                            }

                            auto compute_delta_Gamma = [&](size_t idx, size_t target, size_t low1, size_t low2) -> Real {

                                Real d_dx = Diff::diff_1st(getGammaFunc_X(target, low1, low2), u.x(), idx);
                                
                                Real correction = 0.0L;
                                for(size_t m=0; m<3; ++m) {
                                    Real N_val = conn.getN(m, idx);
                                    if (std::abs(N_val) > 1.0e-10L) { 
                                        Real d_dy = Diff::diff_1st(getGammaFunc_Y(target, low1, low2), u.y(), m);
                                        correction += N_val * d_dy;
                                    }
                                }
                                return d_dx - correction;
                            };

                            Real term_delta_1 = compute_delta_Gamma(k, i, j, l); 
                            Real term_delta_2 = compute_delta_Gamma(l, i, j, k); 
                            Real R_ijkl = term_delta_1 - term_delta_2 + term_comm;

                            sum_R += R_ijkl * factor;
                        }
                    }
                }
                W[i] = sum_R;
            }

            
            Real numerator = innerProduct(W, V);

            return numerator / denominator;
        }
    };

} 
} 
} 

#endif 