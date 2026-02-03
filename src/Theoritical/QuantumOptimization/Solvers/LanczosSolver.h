#ifndef AURELIA_THEORETICAL_QUANTUM_LANCZOS_SOLVER_H
#define AURELIA_THEORETICAL_QUANTUM_LANCZOS_SOLVER_H

#include <vector>
#include <cmath>
#include <functional>
#include <iomanip>

#include "../../../Foundation/LinearAlgebra/DynamicVector.h"

namespace Aurelia {
namespace Theoretical {
namespace QuantumOptimization {
namespace Solvers {

    using namespace Aurelia::Math;
    using HessianOracle = std::function<DynamicVector(const DynamicVector&)>; 

    class LanczosSolver {
    private:
        int max_krylov_dim_;
        
    public:
        LanczosSolver(int max_dim = 20) : max_krylov_dim_(max_dim) {}

        DynamicVector solve(const DynamicVector& g, HessianOracle Hv, Real L) {
            Real beta = g.norm();
            if (beta < 1e-12) return DynamicVector::Zero(g.size());

            int n = g.size();
            int k = std::min((int)n, max_krylov_dim_);

            std::vector<DynamicVector> V;
            V.reserve(k+1);

            std::vector<Real> alphas;
            std::vector<Real> betas;

            DynamicVector v_prev = DynamicVector::Zero(n);
            DynamicVector v_curr = g * (-1.0 / beta); 
            V.push_back(v_curr);
            betas.push_back(beta); 

            for(int j=0; j<k; ++j) {
                DynamicVector w = Hv(v_curr); 
                
                if (j > 0) w = w - v_prev * betas.back(); 
                
                Real alpha = w.dot(v_curr);
                alphas.push_back(alpha);
                
                w = w - v_curr * alpha;
                
                Real beta_next = w.norm();
                
                if (beta_next < 1e-9) { 
                    break; 
                }

                v_prev = v_curr;
                v_curr = w * (1.0 / beta_next);
                V.push_back(v_curr);
                betas.push_back(beta_next);
            }

            DynamicVector y_optimal = solveTridiagonalCubic(alphas, betas, beta, L);

            DynamicVector d = DynamicVector::Zero(n);
            for(size_t i=0; i<y_optimal.size(); ++i) {
                d = d + V[i] * y_optimal[i];
            }
            
            return d;
        }

    private:
        DynamicVector solveTridiagonalCubic(const std::vector<Real>& alpha, 
                                            const std::vector<Real>& beta, 
                                            Real g_norm, Real L) {
            size_t k = alpha.size();
            
            DynamicVector y(k, 0.0);
            
            for(int iter=0; iter<100; ++iter) {
                
                DynamicVector Ty(k, 0.0);
                for(size_t i=0; i<k; ++i) {
                    Ty[i] += alpha[i] * y[i];
                    if (i > 0) Ty[i] += beta[i] * y[i-1];
                    if (i < k-1) Ty[i] += beta[i+1] * y[i+1];
                }

                DynamicVector grad(k, 0.0);
                grad[0] = -g_norm; 
                grad = grad + Ty;
                
                Real y_norm = y.norm();
                Real y_norm_sq = y.squaredNorm();
                
                DynamicVector cubic_grad_term = (std::abs(y_norm) > 1e-12) ? y * (L * 0.5 * y_norm) : DynamicVector::Zero(k);
                grad = grad + cubic_grad_term;
                
                if (grad.norm() < 1e-9) break; 

                std::vector<std::vector<Real>> H_sub(k, std::vector<Real>(k, 0.0));
                
                Real term1 = (L * 0.5 * y_norm);
                Real term2 = (std::abs(y_norm) > 1e-12) ? (L * 0.5 / y_norm) : 0.0;
                
                for(size_t r=0; r<k; ++r) {
                    for(size_t c=0; c<k; ++c) {
                        if (r==c) H_sub[r][c] += alpha[r];
                        else if (r == c+1 || c == r+1) H_sub[r][c] += beta[std::min(r,c)+1]; 
                        
                        if (r==c) H_sub[r][c] += term1;
                        
                        H_sub[r][c] += term2 * y[r] * y[c];
                    }
                }
                
                std::vector<Real> step(k);
                
                auto A = H_sub; 
                auto b_vec = grad;
                
                for(size_t i=0; i<k; ++i) {
                    Real maxEl = std::abs(A[i][i]);
                    size_t maxRow = i;
                    for(size_t r=i+1; r<k; ++r) {
                        if(std::abs(A[r][i]) > maxEl) {
                            maxEl = std::abs(A[r][i]);
                            maxRow = r;
                        }
                    }
                    
                    std::swap(A[i], A[maxRow]);
                    Real tmp_b = b_vec[i]; b_vec[i] = b_vec[maxRow]; b_vec[maxRow] = tmp_b;
                    
                    if(std::abs(A[i][i]) < 1e-12) continue; 
                    
                    for(size_t r=i+1; r<k; ++r) {
                        Real factor = -A[r][i] / A[i][i];
                        for(size_t c=i; c<k; ++c) {
                            A[r][c] += factor * A[i][c];
                        }
                        b_vec[r] += factor * b_vec[i]; 
                    }
                }
                
                DynamicVector newton_step(k, 0.0);
                for(int i=k-1; i>=0; --i) {
                   Real sum = 0.0;
                   for(size_t j=i+1; j<k; ++j) {
                       sum += A[i][j] * newton_step[j];
                   }
                   newton_step[i] = (b_vec[i] - sum) / A[i][i];
                }
                
                y = y - newton_step;
            }
            return y;
        }
    };

}
}
}
}

#endif
