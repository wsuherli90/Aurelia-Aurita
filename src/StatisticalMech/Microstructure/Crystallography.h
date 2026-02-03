#ifndef AURELIA_STATMECH_MICROSTRUCTURE_CRYSTALLOGRAPHY_H
#define AURELIA_STATMECH_MICROSTRUCTURE_CRYSTALLOGRAPHY_H

#include <vector>
#include <cmath>
#include <algorithm>

#include "../../../Foundation/LinearAlgebra/Tensor4.h"
#include "../../../Foundation/LinearAlgebra/Matrix.h"
#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace StatMech {
namespace Microstructure {

    using Real = long double;
    using Tensor4Elastic = Aurelia::Math::Tensor4<Real, 3>; 
    using Tensor4Optic   = Aurelia::Math::Tensor4<Real, 4>; 

    class Crystallography {
    public:
 
        static Tensor4Elastic buildStiffnessTensor(Real E_L, Real E_T, 
                                                   Real G_LT, Real v_LT, Real v_TT) {
            Tensor4Elastic C;
            

            
            Real delta = (1.0L + v_TT) * (1.0L - v_TT - 2.0L * v_LT * v_LT * (E_T / E_L));

            Real C11 = E_T * (1.0L - v_LT * v_LT * (E_T / E_L)) / delta;
            Real C33 = E_L * (1.0L - v_TT * v_TT) / delta;
            Real C12 = E_T * (v_TT + v_LT * v_LT * (E_T / E_L)) / delta;
            Real C13 = E_T * v_LT * (1.0L + v_TT) / delta;
            Real C44 = G_LT;
            Real C66 = (C11 - C12) / 2.0L;

            auto set_sym = [&](size_t i, size_t j, size_t k, size_t l, Real val) {

                C(i,j,k,l) = val; C(j,i,k,l) = val; 
                C(i,j,l,k) = val; C(j,i,l,k) = val;
 
                C(k,l,i,j) = val; C(l,k,i,j) = val; 
                C(k,l,j,i) = val; C(l,k,j,i) = val;
            };


            set_sym(0,0, 0,0, C11); /
            set_sym(1,1, 1,1, C11); 
            set_sym(2,2, 2,2, C33); 
            
            set_sym(0,0, 1,1, C12); 
            set_sym(0,0, 2,2, C13); 
            set_sym(1,1, 2,2, C13); 

            set_sym(1,2, 1,2, C44); 
            set_sym(0,2, 0,2, C44); 
            set_sym(0,1, 0,1, C66); 

            return C;
        }

   
        static Tensor4Optic buildOpticalTensor(Real n_iso, Real chirality) {
            Tensor4Optic Chi;
            
   
            Real eps = n_iso * n_iso;
            Real mu_inv = 1.0L; /

            for (size_t i = 1; i < 4; ++i) {
                Chi(0, i, 0, i) = eps;
                Chi(i, 0, i, 0) = eps; 
            }


            Chi(1, 2, 1, 2) = mu_inv; Chi(2, 1, 2, 1) = mu_inv;
            Chi(2, 3, 2, 3) = mu_inv; Chi(3, 2, 3, 2) = mu_inv;
            Chi(3, 1, 3, 1) = mu_inv; Chi(1, 3, 1, 3) = mu_inv;


            Chi.enforceElectromagneticSymmetry();

            if (std::abs(chirality) > Aurelia::Config::NUMERICAL_EPSILON) {
                Tensor4Optic AxionTerm = Tensor4Optic::LeviCivita();
                Chi = Chi + (AxionTerm * chirality);
            }

            return Chi;
        }
    };

} 
} 
} 

#endif 