

#ifndef AURELIA_GEOMETRY_CARTAN_TORSION_H
#define AURELIA_GEOMETRY_CARTAN_TORSION_H

#include <vector>
#include <cmath>
#include "MetricTensor.h"

namespace Aurelia {
namespace Geometry {
namespace Manifold {

    using Real = long double;
    using Matrix3 = Aurelia::Math::Matrix<Real>;
    using Vector3 = Aurelia::Math::Vector<Real>;

    class CartanTorsion {
    private:

        std::vector<Matrix3> data_; 

    public:
        CartanTorsion() : data_(3, Matrix3(3, 3)) {}


        void compute(MetricTensor& metric_engine, PointTM& u) {
            Real h = 1.0e-4L; 
            Real inv_12h = 1.0L / (12.0L * h);
            Vector3 original_y = u.y();

            for (size_t k = 0; k < 3; ++k) {
                
    
                Vector3 y_p2 = original_y; y_p2[k] += 2.0L * h;
                u.set_y(y_p2); metric_engine.compute(u);
                Matrix3 g_p2 = metric_engine.covariant();

                Vector3 y_p1 = original_y; y_p1[k] += h;
                u.set_y(y_p1); metric_engine.compute(u);
                Matrix3 g_p1 = metric_engine.covariant();

                Vector3 y_m1 = original_y; y_m1[k] -= h;
                u.set_y(y_m1); metric_engine.compute(u);
                Matrix3 g_m1 = metric_engine.covariant();

                Vector3 y_m2 = original_y; y_m2[k] -= 2.0L * h;
                u.set_y(y_m2); metric_engine.compute(u);
                Matrix3 g_m2 = metric_engine.covariant();

                // Optimized: compute derivative in-place to avoid temporary matrices
                // dg/dy = (-g_p2 + 8*g_p1 - 8*g_m1 + g_m2) / (12*h)
                // C_ijk = 0.5 * dg_ij/dy_k
                for (size_t i = 0; i < 3; ++i) {
                    for (size_t j = 0; j < 3; ++j) {
                        Real dg_ij = (-g_p2(i,j) + 8.0L*g_p1(i,j) - 8.0L*g_m1(i,j) + g_m2(i,j)) * inv_12h;
                        data_[k](i, j) = dg_ij * 0.5L;
                    }
                }
            }


            u.set_y(original_y);
            metric_engine.compute(u); 
        }


        Real operator()(size_t i, size_t j, size_t k) const {
            return data_[k](i, j);
        }

        Real norm() const {
            Real sum = 0.0L;
            for (size_t k = 0; k < 3; ++k) {
                for (size_t i = 0; i < 3; ++i) {
                    for (size_t j = 0; j < 3; ++j) {
                        Real val = data_[k](i, j);
                        sum += val * val;
                    }
                }
            }
            return std::sqrt(sum);
        }
    };

} 
} 
} 

#endif 