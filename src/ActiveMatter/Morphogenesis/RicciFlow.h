#ifndef AURELIA_ACTIVEMATTER_MORPHOGENESIS_RICCI_FLOW_H
#define AURELIA_ACTIVEMATTER_MORPHOGENESIS_RICCI_FLOW_H

#include <vector>
#include <cmath>
#include <omp.h>
#include <array>
#include <algorithm>
#include <functional>

#include "../Thermodynamics/ChemicalPotential.h"
#include "../../Geometry/Manifold/MetricTensor.h"
#include "../../Geometry/Connection/ChernConnection.h"
#include "../../Foundation/LinearAlgebra/Vector.h"
#include "../../Foundation/LinearAlgebra/Matrix.h"
#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace ActiveMatter {
namespace Morphogenesis {

    using Real = long double;
    using Vector3 = Aurelia::Math::Vector<Real, 3>;
    using Matrix3 = Aurelia::Math::Matrix<Real, 3, 3>;
    using PointTM = Aurelia::Geometry::Manifold::PointTM;
    using ChemPot = Aurelia::ActiveMatter::Thermodynamics::ChemicalPotential;

    struct MetricFieldSoA {
        size_t nx, ny, nz;
        size_t total_size;
        
        std::vector<Real> g11, g12, g13, g22, g23, g33;
        std::vector<Vector3> directors;

        MetricFieldSoA(size_t x, size_t y, size_t z) 
            : nx(x), ny(y), nz(z), total_size(x * y * z) {
            g11.resize(total_size); g12.resize(total_size); g13.resize(total_size);
            g22.resize(total_size); g23.resize(total_size); g33.resize(total_size);
            directors.resize(total_size);
        }

        inline void set(size_t idx, const Matrix3& m) {
            g11[idx] = m(0,0); g12[idx] = m(0,1); g13[idx] = m(0,2);
            g22[idx] = m(1,1); g23[idx] = m(1,2); g33[idx] = m(2,2);
        }

        inline Matrix3 get(size_t idx) const {
            if (idx >= g11.size()) return Matrix3(0.0L); 
            Matrix3 m;
            m(0,0) = g11[idx]; m(0,1) = g12[idx]; m(0,2) = g13[idx];
            m(1,0) = g12[idx]; m(1,1) = g22[idx]; m(1,2) = g23[idx];
            m(2,0) = g13[idx]; m(2,1) = g23[idx]; m(2,2) = g33[idx];
            return m;
        }

        MetricFieldSoA combine(const MetricFieldSoA& other, Real scale) const {
             MetricFieldSoA result(nx, ny, nz);
             result.directors = this->directors; 
             
             #pragma omp parallel for
             for(size_t i=0; i<total_size; ++i) {
                 result.g11[i] = g11[i] + scale * other.g11[i];
                 result.g12[i] = g12[i] + scale * other.g12[i];
                 result.g13[i] = g13[i] + scale * other.g13[i];
                 result.g22[i] = g22[i] + scale * other.g22[i];
                 result.g23[i] = g23[i] + scale * other.g23[i];
                 result.g33[i] = g33[i] + scale * other.g33[i];
             }
             return result;
        }

        void accumulate(const MetricFieldSoA& other, Real scale) {
             #pragma omp parallel for
             for(size_t i=0; i<total_size; ++i) {
                 g11[i] += scale * other.g11[i];
                 g12[i] += scale * other.g12[i];
                 g13[i] += scale * other.g13[i];
                 g22[i] += scale * other.g22[i];
                 g23[i] += scale * other.g23[i];
                 g33[i] += scale * other.g33[i];
             }
        }
    };

    using ChernConn = Aurelia::Geometry::Connection::ChernConnection;

    class RicciFlow {
    private:
        Real dt_; 
        Real dx_;
        size_t nx_, ny_, nz_;

    public:
        explicit RicciFlow(size_t nx, size_t ny, size_t nz, Real dx = Aurelia::Config::GRID_DX) 
            : dt_(Aurelia::Config::TIME_DT), dx_(dx), nx_(nx), ny_(ny), nz_(nz) {}

        struct Workspace {
            std::vector<std::array<Matrix3, 3>> gamma_field;
            std::vector<Matrix3> g_inv_field;
            std::vector<Real> mu_field;
            std::vector<Matrix3> R_tensor_field;
            std::vector<Vector3> v_field;

            void resize(size_t total) {
                if (gamma_field.size() != total) {
                    gamma_field.resize(total);
                    g_inv_field.resize(total);
                    mu_field.resize(total);
                    R_tensor_field.resize(total);
                    v_field.resize(total);
                }
            }
        };

        template <typename MetricEngineType>
        static Matrix3 computeRicciTensor(PointTM& u, 
                                          ChernConn& conn, 
                                          MetricEngineType& g_engine) {
            Matrix3 Ricci(0.0L);
            const Real h = Aurelia::Config::FINITE_DIFFERENCE_STEP; 
            const Real inv_12h = 1.0L / (12.0L * h);

            auto getGamma = [&](const Vector3& x_pos, size_t k, size_t i, size_t j) -> Real {
                PointTM temp = u; temp.set_x(x_pos);
                g_engine.compute(temp);
                conn.compute(g_engine, temp);
                return conn(k, i, j);
            };

            Vector3 x = u.x();
            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = 0; j < 3; ++j) {
                    Real sum = 0.0L;
                    for (size_t k = 0; k < 3; ++k) {
                        Vector3 p2 = x; p2[k] += 2.0L*h;
                        Vector3 p1 = x; p1[k] += h;
                        Vector3 m1 = x; m1[k] -= h;
                        Vector3 m2 = x; m2[k] -= 2.0L*h;

                        Real val_p2 = getGamma(p2, k, i, j);
                        Real val_p1 = getGamma(p1, k, i, j);
                        Real val_m1 = getGamma(m1, k, i, j);
                        Real val_m2 = getGamma(m2, k, i, j);
                        
                        Real d_k_Gamma = (-val_p2 + 8.0L*val_p1 - 8.0L*val_m1 + val_m2) * inv_12h;

                        Vector3 p2j = x; p2j[j] += 2.0L*h;
                        Vector3 p1j = x; p1j[j] += h;
                        Vector3 m1j = x; m1j[j] -= h;
                        Vector3 m2j = x; m2j[j] -= 2.0L*h;

                        Real val_p2j = getGamma(p2j, k, i, k);
                        Real val_p1j = getGamma(p1j, k, i, k);
                        Real val_m1j = getGamma(m1j, k, i, k);
                        Real val_m2j = getGamma(m2j, k, i, k);

                        Real d_j_Gamma = (-val_p2j + 8.0L*val_p1j - 8.0L*val_m1j + val_m2j) * inv_12h;
                        
                        conn.compute(g_engine, u); 
                        Real term_quad = 0.0L;
                        for(size_t m=0; m<3; ++m) {
                            term_quad += conn(k, k, m) * conn(m, i, j) 
                                       - conn(k, i, m) * conn(m, k, j);
                        }
                        sum += (d_k_Gamma - d_j_Gamma + term_quad);
                    }
                    Ricci(i, j) = sum;
                }
            }
            return Ricci;
        }

        template <typename MetricEngineType>
        Matrix3 evolve(PointTM& u, 
                       MetricEngineType& g_engine, 
                       ChernConn& conn,
                       const ChemPot& chem_pot,
                       const Vector3& v_flow,
                       const Matrix3& Dv_flow) {
            
            g_engine.compute(u);
            conn.compute(g_engine, u);
            Matrix3 R_ij = computeRicciTensor(u, conn, g_engine);
            Matrix3 g_curr = g_engine.covariant();
            
            return g_curr - (R_ij * 2.0L * dt_);
        }

        inline size_t idx(int i, int j, int k) const {
            int ic = std::max(0, std::min((int)nx_ - 1, i));
            int jc = std::max(0, std::min((int)ny_ - 1, j));
            int kc = std::max(0, std::min((int)nz_ - 1, k));
            return (kc * ny_ + jc) * nx_ + ic;
        }


        Real d_k_g_ab(const MetricFieldSoA& field, int i, int j, int k, int a, int b, int deriv_dir) const {
             Real inv_12dx = 1.0L / (12.0L * dx_);
             
             auto getVal = [&](int ii, int jj, int kk) {
                 size_t id = idx(ii, jj, kk);
                 if((a==0 && b==0)) return field.g11[id];
                 if((a==0 && b==1) || (a==1 && b==0)) return field.g12[id];
                 if((a==0 && b==2) || (a==2 && b==0)) return field.g13[id];
                 if((a==1 && b==1)) return field.g22[id];
                 if((a==1 && b==2) || (a==2 && b==1)) return field.g23[id];
                 if((a==2 && b==2)) return field.g33[id];
                 return 0.0L;
             };

             if (deriv_dir == 0) return (-getVal(i+2, j, k) + 8.0L*getVal(i+1, j, k) - 8.0L*getVal(i-1, j, k) + getVal(i-2, j, k)) * inv_12dx;
             if (deriv_dir == 1) return (-getVal(i, j+2, k) + 8.0L*getVal(i, j+1, k) - 8.0L*getVal(i, j-1, k) + getVal(i, j-2, k)) * inv_12dx;
             if (deriv_dir == 2) return (-getVal(i, j, k+2) + 8.0L*getVal(i, j, k+1) - 8.0L*getVal(i, j, k-1) + getVal(i, j, k-2)) * inv_12dx;
             return 0.0L;
        }

        Real computeGammaUpper(const MetricFieldSoA& field, const Matrix3& g_inv, 
                               int i_grid, int j_grid, int k_grid, 
                               int upper_k, int lower_i, int lower_j) const {
            Real sum = 0.0L;
            for(int d=0; d<3; ++d) {
                Real term1 = d_k_g_ab(field, i_grid, j_grid, k_grid, lower_i, d, lower_j);
                Real term2 = d_k_g_ab(field, i_grid, j_grid, k_grid, lower_j, d, lower_i);
                Real term3 = d_k_g_ab(field, i_grid, j_grid, k_grid, lower_i, lower_j, d);
                
                Real gamma_lower = 0.5L * (term1 + term2 - term3);
                sum += g_inv(upper_k, d) * gamma_lower;
            }
            return sum;
        }
        
        Vector3 gradientScalar(const std::vector<Real>& field, int i, int j, int k) const {
            Vector3 grad;
            Real inv_12dx = 1.0L / (12.0L * dx_);
            
            grad[0] = (-field[idx(i+2, j, k)] + 8.0*field[idx(i+1, j, k)] - 8.0*field[idx(i-1, j, k)] + field[idx(i-2, j, k)]) * inv_12dx;
            grad[1] = (-field[idx(i, j+2, k)] + 8.0*field[idx(i, j+1, k)] - 8.0*field[idx(i, j-1, k)] + field[idx(i, j-2, k)]) * inv_12dx;
            grad[2] = (-field[idx(i, j, k+2)] + 8.0*field[idx(i, j, k+1)] - 8.0*field[idx(i, j, k-1)] + field[idx(i, j, k-2)]) * inv_12dx;
            return grad;
        }

        Matrix3 gradientVector(const std::vector<Vector3>& field, int i, int j, int k) const {
            Matrix3 grad;
            Real inv_12dx = 1.0L / (12.0L * dx_);
            for(int c=0; c<3; ++c) {
                grad(c, 0) = (-field[idx(i+2, j, k)][c] + 8.0*field[idx(i+1, j, k)][c] - 8.0*field[idx(i-1, j, k)][c] + field[idx(i-2, j, k)][c]) * inv_12dx;
                grad(c, 1) = (-field[idx(i, j+2, k)][c] + 8.0*field[idx(i, j+1, k)][c] - 8.0*field[idx(i, j-1, k)][c] + field[idx(i, j-2, k)][c]) * inv_12dx;
                grad(c, 2) = (-field[idx(i, j, k+2)][c] + 8.0*field[idx(i, j, k+1)][c] - 8.0*field[idx(i, j, k-1)][c] + field[idx(i, j, k-2)][c]) * inv_12dx;
            }
            return grad;
        }


        void computeFlowRate(const MetricFieldSoA& field,
                             MetricFieldSoA& rate_out, 
                             Workspace& ws,
                             const ChemPot& chem_pot) {
            
            size_t total = field.total_size;
            ws.resize(total); 

            std::vector<std::array<Matrix3, 3>>& gamma_field = ws.gamma_field;
            std::vector<Matrix3>& g_inv_field = ws.g_inv_field;

            #pragma omp parallel for schedule(dynamic)
            for (size_t id = 0; id < total; ++id) {
                 Matrix3 g = field.get(id);
                 Real det = g.determinant();
                 Matrix3 g_inv;
                 if(std::abs(det) > 1e-12) g_inv = g.inverse();
                 else g_inv = Matrix3::identity();
                 g_inv_field[id] = g_inv;
            }

            #pragma omp parallel for schedule(dynamic)
            for(size_t id=0; id<total; ++id) {
                int k = id / (nx_ * ny_);
                int j = (id % (nx_ * ny_)) / nx_;
                int i = id % nx_;

                for(int upper_k=0; upper_k<3; ++upper_k) {
                    for(int a=0; a<3; ++a) {
                        for(int b=0; b<3; ++b) {
                            gamma_field[id][upper_k](a, b) = computeGammaUpper(field, g_inv_field[id], i, j, k, upper_k, a, b);
                        }
                    }
                }
            }


            std::vector<Real>& mu_field = ws.mu_field;
            std::vector<Matrix3>& R_tensor_field = ws.R_tensor_field;

            #pragma omp parallel for schedule(dynamic)
            for (size_t id = 0; id < total; ++id) {
                int k = id / (nx_ * ny_);
                int j = (id % (nx_ * ny_)) / nx_;
                int i = id % nx_;

                Real inv_12dx = 1.0L / (12.0L * dx_);

                Matrix3 Ricci;

                auto d_gamma = [&](int upper, int lower1, int lower2, int deriv_comp) {
                     auto getG = [&](int ii, int jj, int kk) {
                          return gamma_field[idx(ii, jj, kk)][upper](lower1, lower2);
                     };

                     int i1=i,j1=j,k1=k; int di=0,dj=0,dk=0;
                     if(deriv_comp==0) di=1;
                     if(deriv_comp==1) dj=1;
                     if(deriv_comp==2) dk=1;

                     Real val_p2 = getG(i+2*di, j+2*dj, k+2*dk);
                     Real val_p1 = getG(i+di, j+dj, k+dk);
                     Real val_m1 = getG(i-di, j-dj, k-dk);
                     Real val_m2 = getG(i-2*di, j-2*dj, k-2*dk);
                     
                     return (-val_p2 + 8.0*val_p1 - 8.0*val_m1 + val_m2) * inv_12dx;
                };

                for(int a=0; a<3; ++a) {
                    for(int b=0; b<3; ++b) {
                        Real sum = 0.0L;
                        for(int c=0; c<3; ++c) sum += d_gamma(c, a, b, c);
                        for(int c=0; c<3; ++c) sum -= d_gamma(c, a, c, b);
                        for(int c=0; c<3; ++c) {
                            for(int d=0; d<3; ++d) {
                                sum += gamma_field[id][c](c, d) * gamma_field[id][d](a, b);
                            }
                        }
                        for(int c=0; c<3; ++c) {
                            for(int d=0; d<3; ++d) {
                                sum -= gamma_field[id][c](a, d) * gamma_field[id][d](c, b);
                            }
                        }
                        Ricci(a,b) = sum;
                    }
                }

                Matrix3 g_inv = g_inv_field[id];
                Real R = (g_inv * Ricci).trace();
                R_tensor_field[id] = Ricci;
                
                PointTM u({(Real)i*dx_, (Real)j*dx_, (Real)k*dx_}, field.directors[id]);
                mu_field[id] = chem_pot.compute(u, R);
            }


            std::vector<Vector3>& v_field = ws.v_field;
            Real mobility = Aurelia::Config::TISSUE_MOBILITY;

            #pragma omp parallel for schedule(static)
            for (size_t id = 0; id < total; ++id) {
                int k = id / (nx_ * ny_);
                int j = (id % (nx_ * ny_)) / nx_;
                int i = id % nx_;

                Vector3 grad_mu = gradientScalar(mu_field, i, j, k);
                Matrix3 g_inv = g_inv_field[id];
                Vector3 v;
                for(int a=0; a<3; ++a) {
                    Real sum = 0.0L;
                    for(int b=0; b<3; ++b) sum += g_inv(a,b) * grad_mu[b];
                    v[a] = -mobility * sum;
                }
                v_field[id] = v;
            }

            #pragma omp parallel for schedule(dynamic)
            for (size_t id = 0; id < total; ++id) {
                int k = id / (nx_ * ny_);
                int j = (id % (nx_ * ny_)) / nx_;
                int i = id % nx_;

                Vector3 d_mu = gradientScalar(mu_field, i, j, k);
                Matrix3 H_active;
                
                Real inv_12dx2 = 1.0L / (12.0L * dx_ * dx_);
                Real inv_144dx2 = 1.0L / (144.0L * dx_ * dx_);

                auto d2_mu = [&](int ax, int ay) {
                     if (ax == ay) {
                         int i_off=0, j_off=0, k_off=0;
                         if(ax==0) i_off=1; if(ax==1) j_off=1; if(ax==2) k_off=1;

                         Real f0 = mu_field[idx(i,j,k)];
                         Real f_p1 = mu_field[idx(i+i_off, j+j_off, k+k_off)];
                         Real f_m1 = mu_field[idx(i-i_off, j-j_off, k-k_off)];
                         Real f_p2 = mu_field[idx(i+2*i_off, j+2*j_off, k+2*k_off)];
                         Real f_m2 = mu_field[idx(i-2*i_off, j-2*j_off, k-2*k_off)];

                         return (-f_p2 + 16.0*f_p1 - 30.0*f0 + 16.0*f_m1 - f_m2) * inv_12dx2;
                     } else {
                         auto v = [&](int di, int dj, int dk) {
                             return mu_field[idx(i+di, j+dj, k+dk)];
                         };

                         int di=0, dj=0, dk=0;
                         if(ax==0) di=1; if(ax==1) dj=1; if(ax==2) dk=1; 
                         
                         int qi=0, qj=0, qk=0;
                         if(ay==0) qi=1; if(ay==1) qj=1; if(ay==2) qk=1; 
                         
                         Real sum = 0.0L;
                         int offsets[4] = {2, 1, -1, -2};
                         Real weights[4] = {-1.0L, 8.0L, -8.0L, 1.0L};

                         for(int m=0; m<4; ++m) {
                             for(int n=0; n<4; ++n) {
                                 int oi = offsets[m]*di + offsets[n]*qi;
                                 int oj = offsets[m]*dj + offsets[n]*qj;
                                 int ok = offsets[m]*dk + offsets[n]*qk;
                                 
                                 sum += weights[m] * weights[n] * v(oi, oj, ok);
                             }
                         }
                         return sum * inv_144dx2;
                     }
                };

                for(int a=0; a<3; ++a) {
                    for(int b=0; b<3; ++b) {
                        Real partial_sq = d2_mu(a, b);
                        Real gamma_term = 0.0L;
                        for(int c=0; c<3; ++c) gamma_term += gamma_field[id][c](a, b) * d_mu[c];
                        H_active(a,b) = partial_sq - gamma_term;
                    }
                }

                Vector3 v = v_field[id];
                Matrix3 grad_v = gradientVector(v_field, i, j, k);
                Matrix3 g_curr = field.get(id);
                
                Matrix3 Lie;
                for(int a=0; a<3; ++a) {
                    for(int b=0; b<3; ++b) {
                        Real advection = 0.0L;
                        for(int c=0; c<3; ++c) {
                            advection += v[c] * d_k_g_ab(field, i, j, k, a, b, c);
                        }
                        Real strain = 0.0L;
                        for(int c=0; c<3; ++c) {
                            strain += g_curr(c, b) * grad_v(c, a) + g_curr(a, c) * grad_v(c, b);
                        }
                        Lie(a, b) = advection + strain;
                    }
                }

                Matrix3 R_ij = R_tensor_field[id];
                
                Matrix3 rate = (R_ij * -2.0L) + H_active + Lie;
                
                rate_out.set(id, rate);
            }
        }

    public:
        void evolveField(MetricFieldSoA& field, 
                         const ChemPot& chem_pot,
                         const std::function<Vector3(const PointTM&, Real)>&  ) {
            
            MetricFieldSoA k1(nx_, ny_, nz_);
            MetricFieldSoA k2(nx_, ny_, nz_);
            MetricFieldSoA k3(nx_, ny_, nz_);
            MetricFieldSoA k4(nx_, ny_, nz_);
            
            Workspace ws; 
            ws.resize(field.total_size);

            computeFlowRate(field, k1, ws, chem_pot);

            MetricFieldSoA temp1 = field.combine(k1, 0.5L * dt_);
            computeFlowRate(temp1, k2, ws, chem_pot);

            MetricFieldSoA temp2 = field.combine(k2, 0.5L * dt_);
            computeFlowRate(temp2, k3, ws, chem_pot);

            MetricFieldSoA temp3 = field.combine(k3, dt_);
            computeFlowRate(temp3, k4, ws, chem_pot);

            field.accumulate(k1, dt_ / 6.0L);
            field.accumulate(k2, dt_ / 3.0L); 
            field.accumulate(k3, dt_ / 3.0L);
            field.accumulate(k4, dt_ / 6.0L);
        }
    };

} 
} 
} 

#endif