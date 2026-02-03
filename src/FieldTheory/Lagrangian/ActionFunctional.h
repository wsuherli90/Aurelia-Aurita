#ifndef AURELIA_FIELDTHEORY_LAGRANGIAN_ACTION_FUNCTIONAL_H
#define AURELIA_FIELDTHEORY_LAGRANGIAN_ACTION_FUNCTIONAL_H

#include <functional>
#include <vector>
#include <cmath>

#include "../../Geometry/Manifold/MetricTensor.h"
#include "../../Geometry/Manifold/SlitTangentBundle.h"
#include "../../Foundation/Calculus/Integration.h"
#include "../../Foundation/LinearAlgebra/Vector.h"
#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace FieldTheory {
namespace Lagrangian {

    using Real = long double;
    using Vector3 = Aurelia::Math::Vector<Real>;
    using PointTM = Aurelia::Geometry::Manifold::PointTM;
    using MetricTensor = Aurelia::Geometry::Manifold::MetricTensor;
    using LagrangianDensity = std::function<Real(const PointTM&)>;
    using DirectorField = std::function<Vector3(const Vector3&)>;

    class ActionFunctional {
    private:
        LagrangianDensity L_;
        MetricTensor* g_engine_; 
        DirectorField n_field_;  
        Real x_min_, x_max_;
        Real y_min_, y_max_;
        Real z_min_, z_max_;

    public:
        ActionFunctional(LagrangianDensity lagrangian, MetricTensor* metric)
            : L_(lagrangian), g_engine_(metric),
              x_min_(-0.01L), x_max_(0.01L), 
              y_min_(-0.01L), y_max_(0.01L),
              z_min_(-0.01L), z_max_(0.01L) {
            n_field_ = [](const Vector3&) { return Vector3({0.0L, 0.0L, 1.0L}); };
        }
        void setDomain(Real x0, Real x1, Real y0, Real y1, Real z0, Real z1) {
            x_min_ = x0; x_max_ = x1;
            y_min_ = y0; y_max_ = y1;
            z_min_ = z0; z_max_ = z1;
        }

        void setDirectorField(DirectorField n) {
            n_field_ = n;
        }
        Real compute() {
            using Quad = Aurelia::Math::Calculus::GaussLegendre;
            auto volume_integrand = [&](Real x, Real y, Real z) -> Real {
                Vector3 pos({x, y, z});
                Vector3 dir = n_field_(pos);
                if (dir.squaredNorm() < Aurelia::Config::NUMERICAL_EPSILON) {
                    dir = {0.0L, 0.0L, 1.0L}; 
                }
                PointTM u(pos, dir);
                g_engine_->compute(u);
                Real sqrt_det_g = g_engine_->volumeForm();
                Real lagr_val = L_(u);
                return lagr_val * sqrt_det_g;
            };
            auto integ_z = [&](Real x, Real y) {
                return Quad::integrate([&](Real z){ return volume_integrand(x,y,z); }, z_min_, z_max_);
            };

            auto integ_y = [&](Real x) {
                return Quad::integrate([&](Real y){ return integ_z(x,y); }, y_min_, y_max_);
            };

            Real S = Quad::integrate(integ_y, x_min_, x_max_);
            
            return S;
        }
    };

} 
} 
}

#endif 