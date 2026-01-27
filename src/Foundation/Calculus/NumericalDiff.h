#ifndef AURELIA_FOUNDATION_NUMERICALDIFF_H
#define AURELIA_FOUNDATION_NUMERICALDIFF_H

#include <cmath>
#include <functional>
#include <vector>
#include <stdexcept>

#include "../LinearAlgebra/Vector.h"
#include "../../../config/BiophyiscalConstants.h"

namespace Aurelia {
namespace Math {
namespace Calculus {

    using Real = long double;

    class NumericalDiff {
    private:
        static constexpr Real H = Aurelia::Config::NUMERICAL_DIFF_STEP; 

    public:
        static Real diff_1st(const std::function<Real(const Vector<Real>&)>& func, 
                             Vector<Real> point, 
                             size_t dim_index) {
            
            Real original_val = point[dim_index];

            point[dim_index] = original_val + 2.0 * H;
            Real f_p2 = func(point);

            point[dim_index] = original_val + 1.0 * H;
            Real f_p1 = func(point);

            point[dim_index] = original_val - 1.0 * H;
            Real f_m1 = func(point);

            point[dim_index] = original_val - 2.0 * H;
            Real f_m2 = func(point);


            point[dim_index] = original_val;

            return (-f_p2 + 8.0*f_p1 - 8.0*f_m1 + f_m2) / (12.0 * H);
        }

        static Real diff_2nd(const std::function<Real(const Vector<Real>&)>& func, 
                             Vector<Real> point, 
                             size_t dim_index) {
            
            Real original_val = point[dim_index];
            Real f_0 = func(point); 

            point[dim_index] = original_val + 2.0 * H;
            Real f_p2 = func(point);

            point[dim_index] = original_val + 1.0 * H;
            Real f_p1 = func(point);

            point[dim_index] = original_val - 1.0 * H;
            Real f_m1 = func(point);

            point[dim_index] = original_val - 2.0 * H;
            Real f_m2 = func(point);

            point[dim_index] = original_val;

            return (-f_p2 + 16.0*f_p1 - 30.0*f_0 + 16.0*f_m1 - f_m2) / (12.0 * H * H);
        }

        static Real diff_mixed(const std::function<Real(const Vector<Real>&)>& func, 
                               Vector<Real> point, 
                               size_t i, size_t j) {
            
            if (i == j) return diff_2nd(func, point, i);

            Real val_i = point[i];
            Real val_j = point[j];

            point[i] = val_i + H; point[j] = val_j + H;
            Real f_pp = func(point);


            point[i] = val_i + H; point[j] = val_j - H;
            Real f_pm = func(point);

            point[i] = val_i - H; point[j] = val_j + H;
            Real f_mp = func(point);

            point[i] = val_i - H; point[j] = val_j - H;
            Real f_mm = func(point);

            point[i] = val_i; point[j] = val_j;

            return (f_pp - f_pm - f_mp + f_mm) / (4.0 * H * H);
        }
    };

}
} 
} 

#endif 