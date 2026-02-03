#ifndef AURELIA_FIELDTHEORY_LAGRANGIAN_EULER_LAGRANGE_H
#define AURELIA_FIELDTHEORY_LAGRANGIAN_EULER_LAGRANGE_H

#include <vector>
#include <functional>
#include <cmath>

#include "ActionFunctional.h"
#include "../../Foundation/LinearAlgebra/Vector.h"
#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace FieldTheory {
namespace Lagrangian {

    using Real = long double;
    using VectorX = Aurelia::Math::Vector<Real>;

    class EulerLagrange {
    public:
        static VectorX computeFunctionalDerivative(
            ActionFunctional& action,
            const VectorX& current_field,
            std::function<void(const VectorX&)> set_field_func) {
            
            size_t n = current_field.size();
            VectorX gradient(n);
            

            Real h = Aurelia::Config::VARIATIONAL_PERTURBATION;
            Real dx = Aurelia::Config::GRID_DX;
            Real vol_elem = dx * dx * dx;

            set_field_func(current_field);
            Real S_0 = action.compute();

            VectorX perturbed_field = current_field;

            for (size_t i = 0; i < n; ++i) {
                Real original_val = perturbed_field[i];
                perturbed_field[i] += h;
                
                set_field_func(perturbed_field);
                Real S_plus = action.compute();


                gradient[i] = (S_plus - S_0) / (h * vol_elem);

                perturbed_field[i] = original_val;
            }

            set_field_func(current_field);
            return gradient;
        }

        static bool checkStationarity(const VectorX& gradient) {
            Real norm = gradient.norm();

            return norm < Aurelia::Config::OPTIMIZATION_TOLERANCE;
        }
    };

} 
} 
} 

#endif