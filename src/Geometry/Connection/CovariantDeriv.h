#ifndef AURELIA_GEOMETRY_COVARIANT_DERIV_H
#define AURELIA_GEOMETRY_COVARIANT_DERIV_H

#include <vector>
#include <functional>
#include <cmath>

#include "ChernConnection.h"
#include "../../Foundation/LinearAlgebra/Vector.h"
#include "../../Foundation/LinearAlgebra/Matrix.h"
#include "../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace Geometry {
namespace Connection {

    using Real = long double;
    using Vector3 = Aurelia::Math::Vector<Real>;
    using Matrix3 = Aurelia::Math::Matrix<Real>;


    using VectorField = std::function<Vector3(const Aurelia::Geometry::Manifold::PointTM&)>;

    class CovariantDerivative {
    public:
        static Matrix3 compute(const VectorField& field, 
                               Aurelia::Geometry::Manifold::PointTM& u,
                               const ChernConnection& conn) {
            const Real h = Aurelia::Config::FINITE_DIFFERENCE_STEP;
            const Real inv_2h = 1.0L / (2.0L * h);
            
            Matrix3 D_V(3, 3);
            Vector3 V = field(u); 
            std::vector<Vector3> dV_dx(3);
            std::vector<Vector3> dV_dy(3);

            Vector3 x_orig = u.x();
            Vector3 y_orig = u.y();


            for(size_t k=0; k<3; ++k) {
                Vector3 x_p = x_orig; x_p[k] += h; u.set_x(x_p); Vector3 Vp = field(u);
                Vector3 x_m = x_orig; x_m[k] -= h; u.set_x(x_m); Vector3 Vm = field(u);
                dV_dx[k] = (Vp - Vm) * inv_2h;
            }
            u.set_x(x_orig); 
            for(size_t m=0; m<3; ++m) {
                Vector3 y_p = y_orig; y_p[m] += h; u.set_y(y_p); Vector3 Vp = field(u);
                Vector3 y_m = y_orig; y_m[m] -= h; u.set_y(y_m); Vector3 Vm = field(u);
                dV_dy[m] = (Vp - Vm) * inv_2h;
            }
            u.set_y(y_orig); 
            for(size_t i=0; i<3; ++i) {     
                for(size_t k=0; k<3; ++k) { 
                    Real term_horizontal = dV_dx[k][i];
                    for(size_t m=0; m<3; ++m) {
                        term_horizontal -= conn.getN(m, k) * dV_dy[m][i];
                    }
                    Real term_connection = 0.0L;
                    for(size_t j=0; j<3; ++j) {
                        term_connection += conn(i, j, k) * V[j];
                    }

                    D_V(i, k) = term_horizontal + term_connection;
                }
            }

            return D_V;
        }
        static Real divergence(const VectorField& field, 
                               Aurelia::Geometry::Manifold::PointTM& u,
                               const ChernConnection& conn) {
            Matrix3 D_V = compute(field, u, conn);
            return D_V.trace(); 
        }
    };

} 
} 
} 

#endif 