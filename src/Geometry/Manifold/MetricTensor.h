#ifndef AURELIA_GEOMETRY_METRIC_TENSOR_H
#define AURELIA_GEOMETRY_METRIC_TENSOR_H

#include <cmath>
#include <stdexcept>
#include <array>
#include <type_traits>

#include "SlitTangentBundle.h"
#include "../../Foundation/LinearAlgebra/Matrix.h"
#include "../../Foundation/LinearAlgebra/Vector.h"
#include "../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace Geometry {
namespace Manifold {

    using Real = long double;
    using Matrix3 = Aurelia::Math::Matrix<Real, 3, 3>;
    using Vector3 = Aurelia::Math::Vector<Real, 3>;

    template <typename FinslerFunc>
    class MetricTensor {
    private:
        FinslerFunc F_;     
        
        Matrix3 g_;        
        Matrix3 g_inv_;      
        Real det_;           
        bool computed_;


        static constexpr Real COEFF_1 = 16.0L;
        static constexpr Real COEFF_2 = 30.0L;
        static constexpr Real DIVISOR_SQ = 12.0L; 

    public:

        explicit MetricTensor(FinslerFunc finsler_func) 
            : F_(std::move(finsler_func)), computed_(false) {}

        void compute(const PointTM& u) {
            auto EnergyLagrangian = [&](const Vector3& y_vec) -> Real {
                PointTM temp_u(u.x(), y_vec); 
                Real val = F_(temp_u);
                return 0.5L * val * val;
            };

            g_ = Matrix3(0.0L);
            Vector3 current_y = u.y();

            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = i; j < 3; ++j) {
                    Real val;
                    if (i == j) {
                        val = diff_2nd_optimized(EnergyLagrangian, current_y, i);
                    } else {
                        val = diff_mixed_optimized(EnergyLagrangian, current_y, i, j);
                    }
                    g_(i, j) = val;
                    g_(j, i) = val; 
                }
            }

            det_ = g_.determinant();


            if (det_ <= Aurelia::Config::NUMERICAL_EPSILON) {
                throw std::runtime_error("MetricTensor: Singularity detected (det <= 0). "
                                         "Violates strong convexity condition.");
            }

            g_inv_ = g_.inverse();
            computed_ = true;
        }

        const Matrix3& covariant() const {
            #ifdef DEBUG_MATH
            if (!computed_) throw std::runtime_error("MetricTensor: Data not computed.");
            #endif
            return g_;
        }

        const Matrix3& contravariant() const {
            #ifdef DEBUG_MATH
            if (!computed_) throw std::runtime_error("MetricTensor: Data not computed.");
            #endif
            return g_inv_;
        }

        Real volumeForm() const {
            #ifdef DEBUG_MATH
            if (!computed_) throw std::runtime_error("MetricTensor: Data not computed.");
            #endif
            return std::sqrt(det_);
        }


        Real evaluateEnergy(const PointTM& u) const {
            return F_(u);
        }

    private:
        template <typename Func>
        inline Real diff_2nd_optimized(Func&& func, Vector3 p, size_t idx) const {
            Real h = Aurelia::Config::FINITE_DIFFERENCE_STEP;
            Real h_sq = h * h;
            Real inv_12h2 = 1.0L / (DIVISOR_SQ * h_sq);

            Real f0 = func(p);
            
            Real original = p[idx];

            p[idx] = original - h;       Real fm1 = func(p);
            p[idx] = original + h;       Real fp1 = func(p);
            p[idx] = original - 2.0L*h;  Real fm2 = func(p);
            p[idx] = original + 2.0L*h;  Real fp2 = func(p);

            Real term = -fp2 + COEFF_1*fp1 - COEFF_2*f0 + COEFF_1*fm1 - fm2;
            return term * inv_12h2;
        }

        template <typename Func>
        inline Real diff_mixed_optimized(Func&& func, Vector3 p, size_t i, size_t j) const {
            Real h = Aurelia::Config::FINITE_DIFFERENCE_STEP;
            Real inv_4h2 = 0.25L / (h * h); 
            
            Real v_i = p[i];
            Real v_j = p[j];

            p[i] = v_i + h; p[j] = v_j + h; Real f_pp = func(p);
            p[i] = v_i + h; p[j] = v_j - h; Real f_pm = func(p);
            p[i] = v_i - h; p[j] = v_j + h; Real f_mp = func(p);
            p[i] = v_i - h; p[j] = v_j - h; Real f_mm = func(p);

            return (f_pp - f_pm - f_mp + f_mm) * inv_4h2;
        }
    };

} 
} 
} 

#endif