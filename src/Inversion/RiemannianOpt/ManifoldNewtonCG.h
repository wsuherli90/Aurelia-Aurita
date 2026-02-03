#ifndef AURELIA_INVERSION_RIEMANNIAN_MANIFOLD_NEWTON_CG_H
#define AURELIA_INVERSION_RIEMANNIAN_MANIFOLD_NEWTON_CG_H

#include <cmath>
#include <iostream>
#include <functional>
#include <algorithm>

#include "../../../Foundation/LinearAlgebra/Vector.h"
#include "../../../Foundation/LinearAlgebra/Matrix.h"
#include "../../../Geometry/Invariants/GeodesicFlow.h"
#include "../../../Geometry/Manifold/MetricTensor.h"
#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace Inversion {
namespace RiemannianOpt {

    using Real = long double;
    using VectorX = Aurelia::Math::Vector<Real>;
    using Matrix3 = Aurelia::Math::Matrix<Real>;
    using PointTM = Aurelia::Geometry::Manifold::PointTM;
    using MetricEngine = Aurelia::Geometry::Manifold::MetricTensor;
    using Geodesic = Aurelia::Geometry::Invariants::GeodesicFlow;
    using ChernConn = Aurelia::Geometry::Connection::ChernConnection;

    class ManifoldNewtonCG {
    private:
        Real trust_region_radius_;
        size_t max_inner_iter_;

    public:

        ManifoldNewtonCG() 
            : trust_region_radius_(Aurelia::Config::OPTIMIZATION_INITIAL_TRUST_RADIUS),
              max_inner_iter_(Aurelia::Config::OPTIMIZATION_MAX_CG_ITER) {}

     
        Real innerProduct(const VectorX& u, const VectorX& v, const Matrix3& g) const {
            Real sum = 0.0L;
            for(size_t i=0; i<3; ++i) {
                for(size_t j=0; j<3; ++j) {
                    sum += g(i, j) * u[i] * v[j];
                }
            }
            return sum;
        }


        VectorX solveInnerSubproblem(const VectorX& grad, 
                                     std::function<VectorX(const VectorX&)> HessianOp,
                                     const MetricEngine& metric) {

            Matrix3 g = metric.covariant();

            VectorX eta(grad.size(), 0.0L); 
            VectorX r = grad;  
            VectorX d = -r;   
            

            Real r_norm_sq = innerProduct(r, r, g); 

            for(size_t j=0; j<max_inner_iter_; ++j) {

                VectorX Hd = HessianOp(d);
                Real kappa = innerProduct(d, Hd, g);


                if (kappa <= 0.0L) {
                    Real d_norm = std::sqrt(innerProduct(d, d, g));
                    if (d_norm > Aurelia::Config::NUMERICAL_EPSILON)
                        return d * (trust_region_radius_ / d_norm); 
                    return d;
                }

                Real alpha = r_norm_sq / kappa;
                VectorX eta_next = eta + d * alpha;

                Real eta_next_norm = std::sqrt(innerProduct(eta_next, eta_next, g));
                if (eta_next_norm >= trust_region_radius_) {

                    return eta_next * (trust_region_radius_ / eta_next_norm);
                }


                eta = eta_next;
                VectorX r_next = r + Hd * alpha;
                
                Real r_next_norm_sq = innerProduct(r_next, r_next, g);
                if (std::sqrt(r_next_norm_sq) < Aurelia::Config::OPTIMIZATION_TOLERANCE) {
                    break; 
                }

                Real beta = r_next_norm_sq / r_norm_sq;
                d = -r_next + d * beta;
                r = r_next;
                r_norm_sq = r_next_norm_sq;
            }

            return eta;
        }


        PointTM step(PointTM x, 
                     const VectorX& grad, 
                     std::function<VectorX(const VectorX&)> HessianOp,
                     MetricEngine& metric,
                     ChernConn& conn) {
            

            metric.compute(x);


            VectorX eta = solveInnerSubproblem(grad, HessianOp, metric);
            PointTM launch_state(x.x(), eta);
            PointTM x_new = Geodesic::exponentialMap(metric, conn, launch_state);
            return x_new;
        }
        
        void setTrustRegionRadius(Real delta) { trust_region_radius_ = delta; }
    };

} 
} 
} 

#endif 