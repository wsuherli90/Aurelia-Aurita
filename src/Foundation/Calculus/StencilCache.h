#ifndef AURELIA_FOUNDATION_CALCULUS_STENCIL_CACHE_H
#define AURELIA_FOUNDATION_CALCULUS_STENCIL_CACHE_H

#include <array>
#include <cmath>
#include "../LinearAlgebra/Matrix.h"
#include "../LinearAlgebra/Vector.h"
#include "../../Geometry/Manifold/SlitTangentBundle.h"
#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace Math {
namespace Calculus {

    using Real = long double;
    using Matrix3 = Aurelia::Math::Matrix<Real, 3, 3>;
    using Vector3 = Aurelia::Math::Vector<Real, 3>;
    using PointTM = Aurelia::Geometry::Manifold::PointTM;

    /**
     * @brief High-Performance Stencil Cache for 4th-Order Finite Differences.
     * * Stores pre-computed Metric Tensors at required grid points.
     * Eliminates redundant WLC evaluations by >90%.
     * * Pattern (13 evaluations total):
     * - Center (1)
     * - Axis 0: -2h, -h, +h, +2h (4)
     * - Axis 1: -2h, -h, +h, +2h (4)
     * - Axis 2: -2h, -h, +h, +2h (4)
     */
    class MetricStencilCache {
    public:
        enum Direction { DIR_X = 0, DIR_Y = 1, DIR_Z = 2 };
        enum VariationMode { SPATIAL_X, VERTICAL_Y };

    private:
        // Storage: [Direction][OffsetIndex]
        // OffsetIndex map: 0->-2h, 1->-h, 2->+h, 3->+2h
        std::array<std::array<Matrix3, 4>, 3> neighbors_;
        Matrix3 center_;
        Real step_size_;
        bool initialized_;

    public:
        MetricStencilCache() : step_size_(0.0), initialized_(false) {}

        /**
         * @brief Populates the cache. THIS IS THE EXPENSIVE PART.
         * Runs exactly ONCE per voxel/operation.
         */
        template <typename MetricEngine>
        void populate(MetricEngine& engine, const PointTM& u, VariationMode mode) {
            step_size_ = Aurelia::Config::FINITE_DIFFERENCE_STEP;
            Real h = step_size_;

            // 1. Compute Center
            // Ideally pass by copy if engine modifies state, but compute() resets state usually.
            PointTM p = u; 
            engine.compute(p);
            center_ = engine.covariant();

            // 2. Compute Neighbors
            // Loop unrolled logic for 3 axes
            for (int axis = 0; axis < 3; ++axis) {
                // Offsets: -2h, -h, +h, +2h
                Real offsets[4] = { -2.0L*h, -1.0L*h, 1.0L*h, 2.0L*h };

                for (int k = 0; k < 4; ++k) {
                    p = u; // Reset to center
                    Vector3 vec = (mode == SPATIAL_X) ? p.x() : p.y();
                    
                    vec[axis] += offsets[k];
                    
                    if (mode == SPATIAL_X) p.set_x(vec);
                    else p.set_y(vec);

                    // HEAVY COMPUTATION HERE (WLC Integral)
                    engine.compute(p); 
                    neighbors_[axis][k] = engine.covariant();
                }
            }
            initialized_ = true;
        }

        /**
         * @brief Calculates 4th-Order Partial Derivative of the Metric Tensor.
         * Operation: d(g_ab) / d(coord_k)
         * Cost: Cheap arithmetic (stencil lookup).
         */
        Real partialDerivative(size_t component_a, size_t component_b, size_t deriv_axis) const {
            // Stencil Coefficients for 4th Order Central Difference:
            // f'(x) = (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / 12h
            
            const auto& axis_data = neighbors_[deriv_axis];
            
            // Map: 0->-2h, 1->-h, 2->+h, 3->+2h
            Real f_m2 = axis_data[0](component_a, component_b);
            Real f_m1 = axis_data[1](component_a, component_b);
            Real f_p1 = axis_data[2](component_a, component_b);
            Real f_p2 = axis_data[3](component_a, component_b);

            Real num = -f_p2 + 8.0L*f_p1 - 8.0L*f_m1 + f_m2;
            return num / (12.0L * step_size_);
        }

        const Matrix3& getCenter() const { return center_; }
    };

}
}
}

#endif