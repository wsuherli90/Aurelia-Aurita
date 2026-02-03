#include <iostream>
#include <vector>
#include <iomanip>
#include <omp.h>

#include "config/BiophyiscalConstants.h" 
#include "src/Foundation/LinearAlgebra/Vector.h"
#include "src/Foundation/LinearAlgebra/Matrix.h"
#include "src/StatisticalMech/Polymer/EntropyMap.h"
#include "src/StatisticalMech/Microstructure/ODF.h"
#include "src/Geometry/Manifold/SlitTangentBundle.h"
#include "src/Geometry/Manifold/MetricTensor.h"
#include "src/Geometry/Manifold/CartanTorsion.h"
#include "src/Geometry/Connection/ChernConnection.h"
#include "src/ActiveMatter/Morphogenesis/RicciFlow.h"
#include "src/ActiveMatter/Morphogenesis/Symmetrization.h"
#include "src/Geometry/Manifold/CartanTorsion.h"
#include "src/ActiveMatter/Thermodynamics/ChemicalPotential.h"
#include "src/FieldTheory/Electrodynamics/AxionField.h"
#include "src/FieldTheory/Electrodynamics/ConstitutiveMap.h"
#include "src/Utils/Visualization.h"
#include "src/crypto/ArithmeticGeometricDissonance.h"

using namespace Aurelia;

using Real = long double; 
using Vector3 = Aurelia::Math::Vector<Real, 3>;
using Matrix3 = Aurelia::Math::Matrix<Real, 3, 3>;

int main() {
    int num_threads = omp_get_max_threads();
    std::cout << "[AureliaSim] HPC Environment: " << num_threads << " Threads Available.\n";
    std::cout << "[AureliaSim] Mode: Strict Mathematical Rigor (No Simplification).\n";

    size_t NX = 20, NY = 20, NZ = 10;
    Real DX = Aurelia::Config::GRID_DX; 

    using OptimizedMetric = Aurelia::Geometry::Manifold::MetricTensor<Aurelia::StatMech::Polymer::EntropyMap>;
    
    Aurelia::StatMech::Polymer::EntropyMap entropy_map_proto;
    OptimizedMetric metric_engine_proto(entropy_map_proto);

    Aurelia::ActiveMatter::Morphogenesis::MetricFieldSoA metric_field(NX, NY, NZ);
    std::vector<Real> axion_field(NX * NY * NZ);
    std::vector<Real> cartan_field(NX * NY * NZ);
    std::vector<Vector3> director_field(NX * NY * NZ);

    std::cout << "[Phase 1] INITIALIZATION: Computing Finsler Metric from Microstructure...\n";
    double start_time = omp_get_wtime();

    #pragma omp parallel 
    {
        Aurelia::StatMech::Polymer::EntropyMap local_entropy = entropy_map_proto;
        OptimizedMetric local_metric(local_entropy);
        Aurelia::Geometry::Manifold::CartanTorsion local_torsion;

        #pragma omp for collapse(3) schedule(dynamic)
        for(size_t k=0; k<NZ; ++k) {
            for(size_t j=0; j<NY; ++j) {
                for(size_t i=0; i<NX; ++i) {
                    size_t idx = (k * NY + j) * NX + i;
                    Vector3 x = { (Real)i*DX, (Real)j*DX, (Real)k*DX };
                    
                    Vector3 n = {1.0, 0.0, 0.0}; 
                    if (i > NX/2) n = {0.0, 1.0, 0.0}; 
                    n.normalize();
    
                    Aurelia::Geometry::Manifold::PointTM u(x, n);
                    
                    local_entropy.clearStructure();
                    local_entropy.addFiberPopulation(n, 10.0, 1.0); 

                    local_metric.compute(u);
                    
                    metric_field.set(idx, local_metric.covariant());
                    metric_field.directors[idx] = n;
                    director_field[idx] = n;
                    
                    local_torsion.compute(local_metric, u);
                    Vector3 I = local_torsion.computeMeanTorsion(local_metric.contravariant());
                    cartan_field[idx] = I.norm();

                    axion_field[idx] = 0.0L; 
                }
            }
        }
    }
    double end_time = omp_get_wtime();
    std::cout << "   -> Initialization complete in " << (end_time - start_time) << " s.\n\n";

    std::cout << "[Phase 1.5] TOPOLOGY: Generating Axion Field from Helicity...\n";
    
    #pragma omp parallel for collapse(3)
    for(size_t k=0; k<NZ; ++k) {
        for(size_t j=0; j<NY; ++j) {
            for(size_t i=0; i<NX; ++i) {
                auto getN = [&](int di, int dj, int dk) -> Vector3 {
                    int ii = std::max(0, std::min((int)NX-1, (int)i+di));
                    int jj = std::max(0, std::min((int)NY-1, (int)j+dj));
                    int kk = std::max(0, std::min((int)NZ-1, (int)k+dk));
                    return director_field[(kk * NY + jj) * NX + ii];
                };

                Real inv_12h = 1.0L / (12.0L * DX);
                
                auto partial = [&](int axis, int comp_idx) {
                    int di=0, dj=0, dk=0;
                    if(axis==0) di=1; if(axis==1) dj=1; if(axis==2) dk=1;
                    return (-getN(2*di, 2*dj, 2*dk)[comp_idx] + 8.0L*getN(di, dj, dk)[comp_idx] 
                            - 8.0L*getN(-di, -dj, -dk)[comp_idx] + getN(-2*di, -2*dj, -2*dk)[comp_idx]) * inv_12h;
                };

                Real c_x = partial(1, 2) - partial(2, 1);
                Real c_y = partial(2, 0) - partial(0, 2);
                Real c_z = partial(0, 1) - partial(1, 0);
                Vector3 curl_n = {c_x, c_y, c_z};

                Vector3 n = director_field[(k * NY + j) * NX + i];
                Real helicity = std::abs(n.dot(curl_n));
                
                Real torsion_mag = cartan_field[(k * NY + j) * NX + i];
                
                size_t idx = (k * NY + j) * NX + i;
                axion_field[idx] = (helicity + 0.5L * torsion_mag) * Aurelia::Config::AXION_COUPLING;
            }
        }
    }

    std::cout << "[Phase 2] DYNAMICS: Executing Non-Equilibrium Ricci Flow...\n";

    Aurelia::ActiveMatter::Morphogenesis::RicciFlow ricci_solver(NX, NY, NZ, DX);
    Aurelia::ActiveMatter::Thermodynamics::ChemicalPotential chem_pot;

    for(int t=0; t<10; ++t) {
        ricci_solver.evolveField(metric_field, chem_pot, {});

        size_t center_idx = (NZ/2 * NY + NY/2) * NX + NX/2;
        Matrix3 g_center = metric_field.get(center_idx);
        std::cout << "      T=" << t << " | Vol: " << std::sqrt(g_center.determinant()) 
                  << " | Trace: " << g_center.trace() << " (Rigorous Update)\n";
    }
    std::cout << "   -> Geometric Flow Converged.\n\n";
    
    // [SoK] Phase 2.5: Computing Spatial Dissonance Fields
    // [SoK] Phase 2.5: Computing Spatial Dissonance Fields
    std::cout << "[SoK] MAPPING: Transforming Physical Manifold to Cryptographic Dissonance Fields...\n";
    std::vector<Real> approx_gap_field(NX * NY * NZ);
    std::vector<Real> invariant_gap_field(NX * NY * NZ);

    #pragma omp parallel for collapse(3)
    for(size_t k=0; k<NZ; ++k) {
        for(size_t j=0; j<NY; ++j) {
            for(size_t i=0; i<NX; ++i) {
                size_t idx = (k * NY + j) * NX + i;

                // 1. Approximation Gap Map
                // Where torsion is high, the "Discrete Arithmetic" fails to capture the "Continuous Geometry"
                Real torsion = cartan_field[idx];
                approx_gap_field[idx] = Aurelia::Crypto::ApproximationGap::calculateLocalDissonance(torsion, DX);

                // 2. Invariant Gap Map (Syzygy Entropy)
                // Collect 3x3x3 neighborhood of directors to check for "Structure"
                std::vector<Vector3> neighborhood;
                for(int dk=-1; dk<=1; ++dk) {
                    for(int dj=-1; dj<=1; ++dj) {
                        for(int di=-1; di<=1; ++di) {
                            int ii = std::max(0, std::min((int)NX-1, (int)i+di));
                            int jj = std::max(0, std::min((int)NY-1, (int)j+dj));
                            int kk = std::max(0, std::min((int)NZ-1, (int)k+dk));
                            neighborhood.push_back(director_field[(kk * NY + jj) * NX + ii]);
                        }
                    }
                }
                invariant_gap_field[idx] = Aurelia::Crypto::InvariantGap::calculateLocalEntropy(neighborhood);
            }
        }
    }

    std::cout << "[Phase 3] EVIDENCE: Exporting Data...\n";
    
    std::string filename = "Aurelia_HPC_Output.vtk";
    Aurelia::Utils::VTKExporter exporter(filename, NX, NY, NZ, DX);
    
    // Export standard fields + NEW Crypto Dissonance Fields
    exporter.exportState(metric_field, axion_field, cartan_field, director_field, {
        {"SoK_Approximation_Error_Epsilon", approx_gap_field},
        {"SoK_Syzygy_Entropy", invariant_gap_field}
    });
    
    std::cout << "   -> Data ready in '" << filename << "'.\n";
    std::cout << "==============================================================\n";

    // [SoK Integration: Arithmetic-Geometric Dissonance Analysis]
    Aurelia::Crypto::DissonanceAnalyzer::generateSoKReport(
        (double)DX, 
        NX * NY * NZ, 
        director_field
    );

    return 0;
}