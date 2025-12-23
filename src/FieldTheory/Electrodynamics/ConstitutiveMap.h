#ifndef AURELIA_FIELDTHEORY_ELECTRODYNAMICS_CONSTITUTIVE_MAP_H
#define AURELIA_FIELDTHEORY_ELECTRODYNAMICS_CONSTITUTIVE_MAP_H
#include <complex>
#include "../../../Foundation/LinearAlgebra/Matrix.h"
#include "../../../Foundation/LinearAlgebra/Tensor4.h"
#include "../../../Foundation/LinearAlgebra/Vector.h"
#include "../../../Foundation/LinearAlgebra/ComplexMath.h"
#include "../../../config/BiophysicalConstants.h"
#include "SkewonField.h"

namespace Aurelia {
namespace FieldTheory {
namespace Electrodynamics {

    using Complex = Aurelia::Math::ComplexH;
    using Real = long double;
    using Tensor4C = Aurelia::Math::Tensor4<Complex, 4>; 
    using Matrix4C = Aurelia::Math::Matrix<Complex>;
    using Vector3 = Aurelia::Math::Vector<Real>;

    class ConstitutiveMap {
    private:
        Tensor4C chi_principal_; 
        Tensor4C chi_total_;     
        SkewonField skewon_engine_;
        Complex current_alpha_;    
        Vector3 current_velocity_; 

        const Tensor4C levi_civita_;

    public:
        ConstitutiveMap() 
            : levi_civita_(Tensor4C::LeviCivita()), 
              skewon_engine_() {
            current_alpha_ = Complex(0.0L);
            current_velocity_ = Vector3(3, 0.0L);
        }


        explicit ConstitutiveMap(const Tensor4C& material_tensor) 
            : chi_principal_(material_tensor), 
              chi_total_(material_tensor), 
              levi_civita_(Tensor4C::LeviCivita()),
              skewon_engine_() {
            current_alpha_ = Complex(0.0L);
            current_velocity_ = Vector3(3, 0.0L);
        }

        void updateState(Complex alpha, const Vector3& velocity) {
            current_alpha_ = alpha;
            current_velocity_ = velocity;
            recomputeTensor();
        }

        void setAxionField(Complex alpha) {
            updateState(alpha, current_velocity_);
        }

 
        void recomputeTensor() {
            chi_total_ = chi_principal_;
            if (current_alpha_.abs() > Aurelia::Config::NUMERICAL_EPSILON) {
                Tensor4C axion_part = levi_civita_ * current_alpha_;
                chi_total_ = chi_total_ + axion_part;
            }
            if (current_velocity_.norm() > Aurelia::Config::NUMERICAL_EPSILON) {
                Tensor4C skewon_part = skewon_engine_.computeTensor(current_velocity_);
                chi_total_ = chi_total_ + skewon_part;
            }
        }

        Matrix4C apply(const Matrix4C& F) const {
            return chi_total_.doubleContract(F);
        }



        Complex getAxionField() const {
            return current_alpha_;
        }

        const Tensor4C& getTensor() const {
            return chi_total_;
        }
    };

} 
} 
} 

#endif 