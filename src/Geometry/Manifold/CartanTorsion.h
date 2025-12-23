

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
            Vector3 original_y = u.y();

            for (size_t k = 0; k < 3; ++k) {
                
    
                Vector3 y_p2 = original_y; y_p2[k] += 2.0 * h;
                u.set_y(y_p2); metric_engine.compute(u);
                Matrix3 g_p2 = metric_engine.covariant();

                Vector3 y_p1 = original_y; y_p1[k] += h;
                u.set_y(y_p1); metric_engine.compute(u);
                Matrix3 g_p1 = metric_engine.covariant();

                Vector3 y_m1 = original_y; y_m1[k] -= h;
                u.set_y(y_m1); metric_engine.compute(u);
                Matrix3 g_m1 = metric_engine.covariant();

                Vector3 y_m2 = original_y; y_m2[k] -= 2.0 * h;
                u.set_y(y_m2); metric_engine.compute(u);
                Matrix3 g_m2 = metric_engine.covariant();

                Matrix3 dg_dy = (g_m2 * -1.0 + g_p2 + g_p1 * 8.0 - g_m1 * 8.0) * (1.0 / (12.0 * h));


                data_[k] = dg_dy * 0.5L;
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