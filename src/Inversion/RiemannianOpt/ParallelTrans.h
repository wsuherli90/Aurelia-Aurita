#ifndef AURELIA_INVERSION_RIEMANNIAN_PARALLEL_TRANS_H
#define AURELIA_INVERSION_RIEMANNIAN_PARALLEL_TRANS_H

#include "../../../Foundation/LinearAlgebra/Vector.h"
#include "../../../Geometry/Manifold/SlitTangentBundle.h"
#include "../../../Geometry/Manifold/MetricTensor.h" 
#include "../../../Geometry/Connection/ChernConnection.h"

namespace Aurelia {
namespace Inversion {
namespace RiemannianOpt {

    using Real = long double;
    using Vector3 = Aurelia::Math::Vector<Real>;
    using PointTM = Aurelia::Geometry::Manifold::PointTM;
    using ChernConn = Aurelia::Geometry::Connection::ChernConnection;
    using MetricTensor = Aurelia::Geometry::Manifold::MetricTensor;

    class ParallelTransport {
    public:
        static Vector3 transport(const Vector3& V_orig, 
                                 const PointTM& start, 
                                 const PointTM& end, 
                                 MetricTensor& g_engine,
                                 ChernConn& conn) {
            

            Vector3 dx = end.x() - start.x();
            

            Vector3 mid_x = (start.x() + end.x()) * 0.5L;
            Vector3 mid_y = (start.y() + end.y()) * 0.5L; 
            
            PointTM mid_point(mid_x, mid_y);
            

            g_engine.compute(mid_point);
            conn.compute(g_engine, mid_point);

            Vector3 V_new = V_orig;

            for (size_t i = 0; i < 3; ++i) {
                Real correction = 0.0L;
                for (size_t j = 0; j < 3; ++j) {
                    for (size_t k = 0; k < 3; ++k) {

                        correction += conn(i, j, k) * V_orig[j] * dx[k];
                    }
                }
                V_new[i] -= correction;
            }

            return V_new;
        }


        static Vector3 inverseTransport(const Vector3& V_at_end, 
                                        const PointTM& start, 
                                        const PointTM& end, 
                                        MetricTensor& g_engine,
                                        ChernConn& conn) {

            
            Vector3 dx = end.x() - start.x();
            Vector3 mid_x = (start.x() + end.x()) * 0.5L;
            Vector3 mid_y = (start.y() + end.y()) * 0.5L;
            PointTM mid_point(mid_x, mid_y);


            g_engine.compute(mid_point);
            conn.compute(g_engine, mid_point);

            Vector3 V_restored = V_at_end;


            for (size_t i = 0; i < 3; ++i) {
                Real correction = 0.0L;
                for (size_t j = 0; j < 3; ++j) {
                    for (size_t k = 0; k < 3; ++k) {
                        correction += conn(i, j, k) * V_at_end[j] * dx[k];
                    }
                }
                V_restored[i] += correction;
            }
            return V_restored;
        }
    };

} 
} 
} 

#endif 