#ifndef AURELIA_FIELDTHEORY_ELECTRODYNAMICS_SKEWON_FIELD_H
#define AURELIA_FIELDTHEORY_ELECTRODYNAMICS_SKEWON_FIELD_H

#include "../../../Foundation/LinearAlgebra/Matrix.h"
#include "../../../Foundation/LinearAlgebra/Tensor4.h"
#include "../../../Geometry/Manifold/SlitTangentBundle.h"
#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace FieldTheory {
namespace Electrodynamics {

    using Real = long double;
    using Complex = Aurelia::Math::ComplexH;
    using Vector3 = Aurelia::Math::Vector<Real>;
    using Tensor4C = Aurelia::Math::Tensor4<Complex, 4>;

    class SkewonField {
    private:
        Real mixing_strength_;

    public:
        SkewonField() {
            mixing_strength_ = Aurelia::Config::SKEWON_MIXING_STRENGTH; 
        }
        Tensor4C computeTensor(const Vector3& flow_velocity) const {
            Tensor4C S_tensor; 
            
            Real vx = flow_velocity[0];
            Real vy = flow_velocity[1];
            Real vz = flow_velocity[2];
            Complex coupling = Complex(0.0L, -1.0L) * mixing_strength_;
            
            auto add_mixing = [&](size_t i, size_t j, Complex val) {

                size_t k_pair, l_pair;
                if (j == 1) { k_pair = 2; l_pair = 3; } 
                else if (j == 2) { k_pair = 3; l_pair = 1; } 
                else { k_pair = 1; l_pair = 2; } 

                S_tensor(0, i, k_pair, l_pair) += val * coupling;
                S_tensor(0, i, l_pair, k_pair) -= val * coupling;
                S_tensor(k_pair, l_pair, 0, i) -= val * coupling; 
                S_tensor(l_pair, k_pair, 0, i) += val * coupling;
            };


            add_mixing(1, 2,  vz); 
            add_mixing(1, 3, -vy); 
            
            add_mixing(2, 1, -vz); 
            add_mixing(2, 3,  vx); 
            
            add_mixing(3, 1,  vy); 
            add_mixing(3, 2, -vx); 

            return S_tensor;
        }
    };

} 
} 
} 

#endif 