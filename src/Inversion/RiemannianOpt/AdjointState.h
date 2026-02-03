#ifndef AURELIA_INVERSION_RIEMANNIAN_ADJOINT_STATE_H
#define AURELIA_INVERSION_RIEMANNIAN_ADJOINT_STATE_H

#include <functional>
#include "../../../Foundation/LinearAlgebra/Vector.h"
#include "../../../FieldTheory/Lagrangian/ActionFunctional.h"

namespace Aurelia {
namespace Inversion {
namespace RiemannianOpt {

    using Real = long double;
    using VectorX = Aurelia::Math::Vector<Real>;

    class AdjointState {
    public:
        static VectorX computeGradient(
            std::function<void()> forward_solver,
            std::function<VectorX(const VectorX&)> adjoint_solver,
            std::function<VectorX(const VectorX&)> sensitivity_kernel,
            const VectorX& data_misfit_source) {
        
            forward_solver();

            VectorX lambda = adjoint_solver(data_misfit_source);

            VectorX gradient = sensitivity_kernel(lambda);

            return gradient;
        }
    };

} 
} 
} 

#endif 