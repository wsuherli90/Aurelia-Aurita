#ifndef AURELIA_FOUNDATION_EIGENSOLVER_H
#define AURELIA_FOUNDATION_EIGENSOLVER_H

#include <cmath>
#include <limits>
#include <utility>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <numeric> 

#include "Matrix.h"
#include "Vector.h"
#include "../../../config/BiophysicalConstants.h"

namespace Aurelia {
namespace Math {


    template <typename T>
    class SelfAdjointEigensolver {
    private:
        Vector<T> eigenvalues_;
        Matrix<T> eigenvectors_;
        bool computed_;

    public:
        SelfAdjointEigensolver() : computed_(false) {}
        explicit SelfAdjointEigensolver(const Matrix<T>& mat) {
            compute(mat);
        }


        void compute(const Matrix<T>& mat) {
            size_t n = mat.rows();
            if (n != mat.cols()) 
                throw std::invalid_argument("Eigensolver: Matrix must be square.");

            Matrix<T> A = mat;
            eigenvectors_ = Matrix<T>::identity(n);
            eigenvalues_ = Vector<T>(n);

            for (size_t i = 0; i < n; ++i) {
                for (size_t j = i + 1; j < n; ++j) {
                    if (std::abs(A(i, j) - A(j, i)) > T(1.0e-5)) { 
                        throw std::invalid_argument("Eigensolver: Matrix is not symmetric.");
                    }
                }
            }

            int max_sweeps = 100; 
            bool converged = false;

            for (int sweep = 0; sweep < max_sweeps; ++sweep) {
                T max_off_diag = T(0);
                for (size_t i = 0; i < n; ++i) {
                    for (size_t j = i + 1; j < n; ++j) {
                        max_off_diag = std::max(max_off_diag, std::abs(A(i, j)));
                    }
                }

                if (max_off_diag < Aurelia::Config::NUMERICAL_EPSILON) {
                    converged = true;
                    break;
                }


                for (size_t p = 0; p < n; ++p) {
                    for (size_t q = p + 1; q < n; ++q) {
                        T app = A(p, p);
                        T aqq = A(q, q);
                        T apq = A(p, q);
                        if (std::abs(apq) < Aurelia::Config::NUMERICAL_EPSILON) continue;
                        T theta = (aqq - app) / (T(2.0) * apq);
                        T t;
                        if (theta >= 0)
                            t = T(1.0) / (theta + std::hypot(T(1.0), theta));
                        else
                            t = T(-1.0) / (-theta + std::hypot(T(1.0), theta));

                        T c = T(1.0) / std::hypot(T(1.0), t);
                        T s = t * c;
                        T tau = s / (T(1.0) + c);
                        A(p, p) -= t * apq;
                        A(q, q) += t * apq;
                        A(p, q) = T(0);
                        A(q, p) = T(0); 

                        for (size_t r = 0; r < n; ++r) {
                            if (r != p && r != q) {
                                T arp = A(r, p);
                                T arq = A(r, q);
                                A(r, p) = c * arp - s * arq;
                                A(p, r) = A(r, p); 
                                A(r, q) = s * arp + c * arq;
                                A(q, r) = A(r, q); 
                            }
                        }
                        for (size_t r = 0; r < n; ++r) {
                            T vrp = eigenvectors_(r, p);
                            T vrq = eigenvectors_(r, q);
                            eigenvectors_(r, p) = c * vrp - s * vrq;
                            eigenvectors_(r, q) = s * vrp + c * vrq;
                        }
                    }
                }
            }

            if (!converged) {
                throw std::runtime_error("Eigensolver: Failed to converge within max sweeps.");
            }
            for (size_t i = 0; i < n; ++i) {
                eigenvalues_[i] = A(i, i);
            }

            sortEigenpairs();
            computed_ = true;
        }

        Vector<T> eigenvalues() const {
            if (!computed_) throw std::runtime_error("Eigensolver: Not computed yet.");
            return eigenvalues_;
        }

        Matrix<T> eigenvectors() const {
            if (!computed_) throw std::runtime_error("Eigensolver: Not computed yet.");
            return eigenvectors_;
        }

    private:

        void sortEigenpairs() {
            size_t n = eigenvalues_.size();
            std::vector<size_t> idx(n);
            std::iota(idx.begin(), idx.end(), 0);
            std::stable_sort(idx.begin(), idx.end(), 
                [this](size_t i1, size_t i2) {
                    return eigenvalues_[i1] > eigenvalues_[i2]; 
                });

            Vector<T> sorted_vals(n);
            Matrix<T> sorted_vecs(n, n);

            for (size_t i = 0; i < n; ++i) {
                sorted_vals[i] = eigenvalues_[idx[i]];
                for (size_t r = 0; r < n; ++r) {
                    sorted_vecs(r, i) = eigenvectors_(r, idx[i]);
                }
            }

            eigenvalues_ = sorted_vals;
            eigenvectors_ = sorted_vecs;
        }
    };

    template <typename T>
    Matrix<T> MatrixExponential(const Matrix<T>& A) {
        SelfAdjointEigensolver<T> solver(A);
        Vector<T> vals = solver.eigenvalues();
        Matrix<T> vecs = solver.eigenvectors();

        Matrix<T> D = Matrix<T>::identity(A.rows());
        for (size_t i = 0; i < vals.size(); ++i) {
            D(i, i) = std::exp(vals[i]);
        }

        return vecs * D * vecs.transpose();
    }


    template <typename T>
    Matrix<T> MatrixLogarithm(const Matrix<T>& A) {
        SelfAdjointEigensolver<T> solver(A);
        Vector<T> vals = solver.eigenvalues();
        Matrix<T> vecs = solver.eigenvectors();

        Matrix<T> D = Matrix<T>::identity(A.rows());
        for (size_t i = 0; i < vals.size(); ++i) {
            if (vals[i] <= Aurelia::Config::NUMERICAL_EPSILON)
                throw std::runtime_error("MatrixLogarithm: Matrix is not Positive Definite.");
            D(i, i) = std::log(vals[i]);
        }

        return vecs * D * vecs.transpose();
    }

    template <typename T>
    Matrix<T> MatrixSqrt(const Matrix<T>& A) {
        SelfAdjointEigensolver<T> solver(A);
        Vector<T> vals = solver.eigenvalues();
        Matrix<T> vecs = solver.eigenvectors();

        Matrix<T> D = Matrix<T>::identity(A.rows());
        for (size_t i = 0; i < vals.size(); ++i) {
            if (vals[i] < T(0)) throw std::runtime_error("MatrixSqrt: Negative Eigenvalue.");
            D(i, i) = std::sqrt(vals[i]);
        }

        return vecs * D * vecs.transpose();
    }
    
 
    template <typename T>
    Matrix<T> MatrixInvSqrt(const Matrix<T>& A) {
        SelfAdjointEigensolver<T> solver(A);
        Vector<T> vals = solver.eigenvalues();
        Matrix<T> vecs = solver.eigenvectors();

        Matrix<T> D = Matrix<T>::identity(A.rows());
        for (size_t i = 0; i < vals.size(); ++i) {
            if (vals[i] <= Aurelia::Config::NUMERICAL_EPSILON) 
                throw std::runtime_error("MatrixInvSqrt: Eigenvalue too close to zero.");
            D(i, i) = T(1.0) / std::sqrt(vals[i]);
        }

        return vecs * D * vecs.transpose();
    }

} 
} 

#endif 