#ifndef AURELIA_FIELDTHEORY_ELECTRODYNAMICS_AXION_FIELD_H
#define AURELIA_FIELDTHEORY_ELECTRODYNAMICS_AXION_FIELD_H

#include <cmath>
#include <functional>

#include "../../../Foundation/LinearAlgebra/ComplexMath.h"
#include "../../../Foundation/LinearAlgebra/Vector.h"
#include "../../../Foundation/Calculus/NumericalDiff.h" 
#include "../../../Geometry/Manifold/SlitTangentBundle.h"
#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace FieldTheory {
namespace Electrodynamics {

    using Real = long double;
    using Complex = Aurelia::Math::ComplexH;
    using Vector3 = Aurelia::Math::Vector<Real>;
    using PointTM = Aurelia::Geometry::Manifold::PointTM;

    class AxionField {
    private:
        Real coupling_strength_;
        std::function<Real(const PointTM&)> chirality_distribution_;

    public:
        AxionField() 
            : coupling_strength_(Aurelia::Config::AXION_COUPLING) { 
            chirality_distribution_ = [](const PointTM&) { return 1.0L; };
        }

        void setDistribution(std::function<Real(const PointTM&)> func) {
            chirality_distribution_ = func;
        }


        Complex evaluate(const PointTM& u) const {
            Real density = chirality_distribution_(u);
            Real real_part = coupling_strength_ * density;
            Real imag_part = coupling_strength_ * density * 0.01L; 

            return Complex(real_part, imag_part);
        }

        Vector3 gradient(const PointTM& u) const {
            using Diff = Aurelia::Math::Calculus::NumericalDiff;
            
            Vector3 grad(3);
            

            auto func_wrapper = [&](const Vector3& x_vec) -> Real {
                PointTM temp = u;
                temp.set_x(x_vec);
                return evaluate(temp).real(); 
            };


            Vector3 x_curr = u.x();
            for(size_t i=0; i<3; ++i) {
                grad[i] = Diff::diff_1st(func_wrapper, x_curr, i);
            }
            return grad;
        }

        Real getCoupling() const { return coupling_strength_; }
    };

} 
} 
} 

#endif 