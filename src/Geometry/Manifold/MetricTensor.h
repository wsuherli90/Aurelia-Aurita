

#ifndef AURELIA_GEOMETRY_METRIC_TENSOR_H
#define AURELIA_GEOMETRY_METRIC_TENSOR_H

#include <functional>
#include <memory>
#include <cmath>
#include <stdexcept>

#include "SlitTangentBundle.h"
#include "../../Foundation/LinearAlgebra/Matrix.h"
#include "../../Foundation/Calculus/NumericalDiff.h"
#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace Geometry {
namespace Manifold {

    using Real = long double;
    using Matrix3 = Aurelia::Math::Matrix<Real>;
    using DiffEngine = Aurelia::Math::Calculus::NumericalDiff;
    using FinslerFunction = std::function<Real(const PointTM&)>;

    class MetricTensor {
    private:
        FinslerFunction F_; 
        
        Matrix3 g_;        
        Matrix3 g_inv_;    
        Real det_;          

        bool computed_;

    public:
        explicit MetricTensor(FinslerFunction finsler_func) 
            : F_(std::move(finsler_func)), computed_(false) {}

        void compute(const PointTM& u) {
            auto EnergyLagrangian = [&](const Aurelia::Math::Vector<Real>& y_vec) -> Real {
                PointTM temp_u(u.x(), y_vec); 
                Real val = F_(temp_u);
                return 0.5L * val * val;
            };

            g_ = Matrix3(3, 3);
            Aurelia::Math::Vector<Real> current_y = u.y();
            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = i; j < 3; ++j) {
                    Real val;
                    if (i == j) {
                        val = DiffEngine::diff_2nd(EnergyLagrangian, current_y, i);
                    } else {
                        val = DiffEngine::diff_mixed(EnergyLagrangian, current_y, i, j);
                    }
                    g_(i, j) = val;
                    g_(j, i) = val; // Enforce Symmetry
                }
            }

            det_ = g_.determinant();

            if (det_ <= Aurelia::Config::NUMERICAL_EPSILON) {
                throw std::runtime_error("MetricTensor: Singularity detected (det <= 0). "
                                         "Material instability or zero-section violation.");
            }

            g_inv_ = g_.inverse();
            
            computed_ = true;
        }


        const Matrix3& covariant() const {
            if (!computed_) throw std::runtime_error("MetricTensor: Data not computed.");
            return g_;
        }

        Real operator()(size_t i, size_t j) const {
             if (!computed_) throw std::runtime_error("MetricTensor: Data not computed.");
             return g_(i, j);
        }

        const Matrix3& contravariant() const {
            if (!computed_) throw std::runtime_error("MetricTensor: Data not computed.");
            return g_inv_;
        }

        Real volumeForm() const {
            if (!computed_) throw std::runtime_error("MetricTensor: Data not computed.");
            return std::sqrt(det_);
        }


        Real evaluateEnergy(const PointTM& u) const {
            return F_(u);
        }
    };

} 
} 
} 

#endif 