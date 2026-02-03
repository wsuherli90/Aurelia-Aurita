#ifndef AURELIA_FIELDTHEORY_ELECTRODYNAMICS_PREMETRIC_MAXWELL_H
#define AURELIA_FIELDTHEORY_ELECTRODYNAMICS_PREMETRIC_MAXWELL_H

#include <vector>
#include <functional>
#include <complex>

#include "ConstitutiveMap.h"
#include "../../Geometry/Connection/CovariantDeriv.h"
#include "../../Foundation/LinearAlgebra/Matrix.h"
#include "../../Foundation/LinearAlgebra/Vector.h"
#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace FieldTheory {
namespace Electrodynamics {

    using Complex = Aurelia::Math::ComplexH;
    using Real = long double;
    using Vector3 = Aurelia::Math::Vector<Real>;
    using Vector4C = std::vector<Complex>; 
    using Matrix4C = Aurelia::Math::Matrix<Complex>; 
    using PointTM = Aurelia::Geometry::Manifold::PointTM;

    class PreMetricMaxwell {
    public:
        static Matrix4C computeFieldStrength(const std::function<Vector4C(const PointTM&)>& A_field, 
                                             PointTM& u, 
                                             Real omega) {
            
            Matrix4C F(4, 4, Complex(0.0));
            Real h = Aurelia::Config::FINITE_DIFFERENCE_STEP; 
            Vector3 x_orig = u.x();
            Vector4C A_center = A_field(u);
            std::vector<Vector4C> grad_A(3, Vector4C(4)); 

            for(size_t i=0; i<3; ++i) { 
                Vector3 xp = x_orig; xp[i] += h; u.set_x(xp); Vector4C Ap = A_field(u);
                Vector3 xm = x_orig; xm[i] -= h; u.set_x(xm); Vector4C Am = A_field(u);
                
                for(size_t nu=0; nu<4; ++nu) {
                    grad_A[i][nu] = (Ap[nu] - Am[nu]) / (2.0L * h);
                }
            }
            u.set_x(x_orig); 
            Complex i_omega = Complex(0.0, -1.0) * omega; 

            for(size_t mu=0; mu<4; ++mu) {
                for(size_t nu=0; nu<4; ++nu) {
                    if (mu == nu) continue;
                    Complex d_mu_A_nu;
                    if (mu == 0) d_mu_A_nu = i_omega * A_center[nu];
                    else d_mu_A_nu = grad_A[mu-1][nu];

                    Complex d_nu_A_mu;
                    if (nu == 0) d_nu_A_mu = i_omega * A_center[mu];
                    else d_nu_A_mu = grad_A[nu-1][mu];

                    F(mu, nu) = d_mu_A_nu - d_nu_A_mu;
                }
            }
            return F; 
        }


        static Vector4C computeMaxwellResidual(const std::function<Matrix4C(const PointTM&)>& H_field, 
                                               const Vector4C& J,
                                               PointTM& u,
                                               Real omega) {
            
            Vector4C residual(4);
            Real h = Aurelia::Config::FINITE_DIFFERENCE_STEP; 
            Vector3 x_orig = u.x();
            Matrix4C H_center = H_field(u);
            Complex i_omega = Complex(0.0, -1.0) * omega;

            for (size_t v = 0; v < 4; ++v) {
                Complex div_H = Complex(0.0);

                div_H += i_omega * H_center(0, v);

                for(size_t i=0; i<3; ++i) {
                    Vector3 xp = x_orig; xp[i] += h; u.set_x(xp); Matrix4C Hp = H_field(u);
                    Vector3 xm = x_orig; xm[i] -= h; u.set_x(xm); Matrix4C Hm = H_field(u);
                    
                    Complex dH_dx = (Hp(i+1, v) - Hm(i+1, v)) / (2.0L * h);
                    div_H += dH_dx;
                }
                
                residual[v] = div_H - J[v];
            }
            u.set_x(x_orig);
            return residual;
        }
    };

} 
} 
} 

#endif