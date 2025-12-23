#ifndef AURELIA_GEOMETRY_GEODESIC_FLOW_H
#define AURELIA_GEOMETRY_GEODESIC_FLOW_H

#include <vector>
#include "../Manifold/SlitTangentBundle.h"
#include "../Connection/ChernConnection.h"
#include "../../config/BiophyiscalConstants.h"

namespace Aurelia {
namespace Geometry {
namespace Invariants {

    using Real = long double;
    using Vector3 = Aurelia::Math::Vector<Real>;
    using namespace Aurelia::Geometry::Manifold;
    using namespace Aurelia::Geometry::Connection;

    class GeodesicFlow {
    private:
        struct State {
            Vector3 x; 
            Vector3 y; 
        };

    public:
        static PointTM exponentialMap(MetricTensor& metric, 
                                      ChernConnection& conn, 
                                      PointTM start_u) {
            
            State current = { start_u.x(), start_u.y() };
            PointTM temp_u = start_u; 

            Real t = 0.0L;
            Real t_end = Aurelia::Config::GEODESIC_TOTAL_TIME;
            Real dt = Aurelia::Config::GEODESIC_TIME_STEP;


            auto get_deriv = [&](const State& s) -> State {
                temp_u.set_x(s.x);
                temp_u.set_y(s.y); 

                metric.compute(temp_u);
                conn.compute(metric, temp_u);

                Vector3 dxdt = s.y;
                Vector3 dydt(3, 0.0L);


                for(size_t i=0; i<3; ++i) {
                    for(size_t j=0; j<3; ++j) {
                        for(size_t k=0; k<3; ++k) {
                            dydt[i] -= conn(i, j, k) * s.y[j] * s.y[k];
                        }
                    }
                }
                return {dxdt, dydt};
            };


            while (t < t_end) {
                State k1 = get_deriv(current);

                State s2 = { current.x + k1.x * (0.5L * dt), current.y + k1.y * (0.5L * dt) };
                State k2 = get_deriv(s2);

                State s3 = { current.x + k2.x * (0.5L * dt), current.y + k2.y * (0.5L * dt) };
                State k3 = get_deriv(s3);

                State s4 = { current.x + k3.x * dt, current.y + k3.y * dt };
                State k4 = get_deriv(s4);

                // Final Update
                current.x = current.x + (k1.x + k2.x*2.0L + k3.x*2.0L + k4.x) * (dt / 6.0L);
                current.y = current.y + (k1.y + k2.y*2.0L + k3.y*2.0L + k4.y) * (dt / 6.0L);

                t += dt;
            }

            return PointTM(current.x, current.y);
        }

        static Vector3 parallelTransport(MetricTensor& metric,
                                         ChernConnection& conn,
                                         PointTM start_u,
                                         Vector3 V_initial) {
            
            struct CoupledState {
                Vector3 x;
                Vector3 y; 
                Vector3 V; 
            };

            CoupledState current = { start_u.x(), start_u.y(), V_initial };
            PointTM temp_u = start_u;

            Real t = 0.0L;
            Real dt = Aurelia::Config::GEODESIC_TIME_STEP; 
            Real t_end = Aurelia::Config::GEODESIC_TOTAL_TIME;

            auto get_deriv = [&](const CoupledState& s) -> CoupledState {
                temp_u.set_x(s.x);
                temp_u.set_y(s.y);
                metric.compute(temp_u);
                conn.compute(metric, temp_u);

                Vector3 dxdt = s.y;
                Vector3 dydt(3, 0.0L); 
                Vector3 dVdt(3, 0.0L); 

                for(size_t i=0; i<3; ++i) {
                    for(size_t j=0; j<3; ++j) {
                        for(size_t k=0; k<3; ++k) {
                            Real Gamma = conn(i, j, k);
                            dydt[i] -= Gamma * s.y[j] * s.y[k];
                            dVdt[i] -= Gamma * s.V[j] * s.y[k]; 
                        }
                    }
                }
                return {dxdt, dydt, dVdt};
            };

            while (t < t_end) {
                CoupledState deriv = get_deriv(current);
                current.x = current.x + deriv.x * dt;
                current.y = current.y + deriv.y * dt;
                current.V = current.V + deriv.V * dt;
                t += dt;
            }

            return current.V;
        }
    };

} /
} 
} 

#endif 