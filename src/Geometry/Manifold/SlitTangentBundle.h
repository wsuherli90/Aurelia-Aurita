
#ifndef AURELIA_GEOMETRY_SLIT_TANGENT_BUNDLE_H
#define AURELIA_GEOMETRY_SLIT_TANGENT_BUNDLE_H

#include <stdexcept>
#include <cmath>
#include <iostream>

#include "../../Foundation/LinearAlgebra/Vector.h"
#include "../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace Geometry {
namespace Manifold {

    using Real = long double;

    using Vector3 = Aurelia::Math::Vector<Real, 3>;

    class PointTM {
    private:
        Vector3 x_; 
        Vector3 y_; 

    public:
        PointTM() : x_(3, 0.0), y_({1.0, 0.0, 0.0}) {}
        PointTM(const Vector3& position, const Vector3& direction) 
            : x_(position), y_(direction) {
            validateSlitCondition();
        }
        void validateSlitCondition() const {
            Real sqNorm = y_.squaredNorm();
            Real eps = Aurelia::Config::NUMERICAL_EPSILON;
            
            if (sqNorm < eps * eps) {
                throw std::runtime_error("Finsler Geometry Error: 'y' vector collapsed to Zero Section. "
                                         "Metric is undefined at y=0.");
            }
        }

        const Vector3& x() const { return x_; }
        const Vector3& y() const { return y_; }

        void set_y(const Vector3& new_y) {
            y_ = new_y;
            validateSlitCondition(); 
        }
        
        void set_x(const Vector3& new_x) {
            x_ = new_x;
        }


        friend std::ostream& operator<<(std::ostream& os, const PointTM& p) {
            os << "State(x=" << p.x_ << ", y=" << p.y_ << ")";
            return os;
        }
    };

}
} 
} 

#endif 