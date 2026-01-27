#ifndef AURELIA_GEOMETRY_CARTAN_TORSION_H
#define AURELIA_GEOMETRY_CARTAN_TORSION_H

#include "../../Foundation/LinearAlgebra/Tensor4.h"
#include "../../Foundation/Calculus/StencilCache.h"
#include "MetricTensor.h"
#include <array>

namespace Aurelia {
namespace Geometry {
namespace Manifold {

    using Real = long double;

    class CartanTorsion {
    private:
        std::array<Real, 27> C_;

    public:
        CartanTorsion() {
            C_.fill(0.0L);
        }

        template <typename MetricEngine>
        void compute(MetricEngine& g_engine, const PointTM& u) {
            
            Aurelia::Math::Calculus::MetricStencilCache stencil;
            stencil.populate(g_engine, u, Aurelia::Math::Calculus::MetricStencilCache::VERTICAL_Y);


            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = i; j < 3; ++j) { 
                    for (size_t k = 0; k < 3; ++k) {
                        
                        Real d_gij_dyk = stencil.partialDerivative(i, j, k);
                        
                        Real val = 0.5L * d_gij_dyk;
                        
                        auto set = [&](size_t a, size_t b, size_t c, Real v) {
                             C_[a*9 + b*3 + c] = v;
                        };
                        
                        set(i, j, k, val); set(i, k, j, val); 
                        set(j, i, k, val); set(j, k, i, val);
                        set(k, i, j, val); set(k, j, i, val);
                    }
                }
            }
        }

        Real operator()(size_t i, size_t j, size_t k) const {
            return C_[i*9 + j*3 + k];
        }

        Vector3 computeMeanTorsion(const Matrix3& g_inv) const {
            Vector3 I = {0.0L, 0.0L, 0.0L};
            for(size_t k=0; k<3; ++k) {
                for(size_t i=0; i<3; ++i) {
                    for(size_t j=0; j<3; ++j) {
                        I[k] += g_inv(i, j) * (*this)(i, j, k);
                    }
                }
            }
            return I;
        }
    };

}
}
}

#endif