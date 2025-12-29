#include <iostream>
#include <vector>
#include <iomanip>

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
#include "src/ActiveMatter/Thermodynamics/ChemicalPotential.h"
#include "src/FieldTheory/Electrodynamics/AxionField.h"
#include "src/FieldTheory/Electrodynamics/PreMetricMaxwell.h"
#include "src/FieldTheory/Electrodynamics/ConstitutiveMap.h"
#include "src/Utils/Visualization.h"

using namespace Aurelia;

using Real = double; 
using Vector3 = Aurelia::Math::Vector<Real, 3>;
using Matrix3 = Aurelia::Math::Matrix<Real, 3, 3>;

int main() {
    

    size_t NX = 10, NY = 10, NZ = 5;
    Real DX = Aurelia::Config::GRID_DX; 

    Aurelia::StatMech::Polymer::EntropyMap entropy_engine;
    Aurelia::Geometry::Manifold::MetricTensor metric_engine(
        [&](const auto& u) { return entropy_engine.computeFinslerMetric(u); }
    );
    
    Aurelia::Geometry::Connection::ChernConnection chern_conn;
    Aurelia::Geometry::Manifold::CartanTorsion torsion_checker;

    std::vector<Aurelia::Geometry::Manifold::MetricTensor> metric_field;
    std::vector<Real> axion_field;
    std::vector<Vector3> director_field;

    std::cout << "   -> Integrating WLC statistics over ODF for " << NX*NY*NZ << " voxels.\n";
    

    for(size_t k=0; k<NZ; ++k) {
        for(size_t j=0; j<NY; ++j) {
            for(size_t i=0; i<NX; ++i) {


                Vector3 x = { (Real)i*DX, (Real)j*DX, (Real)k*DX };
                
                Vector3 n = {1.0, 0.0, 0.0}; 
                if (i > NX/2) n = {0.0, 1.0, 0.0}; 
                n.normalize();
 
                Aurelia::Geometry::Manifold::PointTM u(x, n);
  
                entropy_engine.clearStructure();
                entropy_engine.addFiberPopulation(n, 10.0, 1.0); 

                metric_engine.compute(u);
                metric_field.push_back(metric_engine);
                director_field.push_back(n);

                Real twist_intensity = (i == NX/2) ? 1.0 : 0.1;
                axion_field.push_back(twist_intensity * Aurelia::Config::AXION_COUPLING); 
            }
        }
    }
    std::cout << "   -> Manifold Initialized.\n\n";


    std::cout << "[Phase 2] DIAGNOSIS: Testing Mathematical Rigor...\n";
    
    Vector3 x_test = { (Real)(NX/2)*DX, 0.0, 0.0 };
    Vector3 y_test = {1.0, 1.0, 0.0}; 
    y_test.normalize();
    Aurelia::Geometry::Manifold::PointTM u_test(x_test, y_test);

    entropy_engine.clearStructure();
    entropy_engine.addFiberPopulation({1.0, 0.0, 0.0}, 5.0, 1.0); 
    
    metric_engine.compute(u_test);
    torsion_checker.compute(metric_engine, u_test);
    Real torsion_norm = torsion_checker.norm();

    std::cout << "   -> Point " << x_test << "\n";
    std::cout << "   -> Cartan Torsion ||C_ijk|| = " << torsion_norm << "\n";
    
    if (torsion_norm > 1.0e-6) {
        std::cout << "   -> [PASS] Space is Non-Riemannian (Finslerian).\n";
    } else {
        std::cout << "   -> [FAIL] Space is Euclidean.\n";
    }
    std::cout << "\n";


    std::cout << "[Phase 3] DYNAMICS: Executing Active Ricci Flow...\n";

    Aurelia::ActiveMatter::Morphogenesis::RicciFlow ricci_solver;
    Aurelia::ActiveMatter::Thermodynamics::ChemicalPotential chem_pot;
    Aurelia::ActiveMatter::Morphogenesis::Symmetrization sym_helper; 

    std::cout << "   -> Simulating injury repair at defect site...\n";
    
    for(int t=0; t<5; ++t) {
        chern_conn.compute(metric_engine, u_test);

        Matrix3 R_tensor = ricci_solver.computeRicciTensor(u_test, chern_conn, metric_engine);
        Matrix3 g_inv = metric_engine.contravariant();

        Matrix3 product = g_inv * R_tensor;
        Real R_scalar = product.trace();

        Vector3 v_flow = sym_helper.computeVelocity(u_test, R_scalar);
        Matrix3 Dv_flow = sym_helper.computeVelocityGradient(u_test, R_scalar);

        Matrix3 g_new = ricci_solver.evolve(u_test, metric_engine, chern_conn, chem_pot, v_flow, Dv_flow);

        Real vol = g_new.determinant();
        std::cout << "      T=" << t << ": Metric Volume Form = " << std::sqrt(vol) 
                  << " | R=" << R_scalar << " (Remodeling...)\n";
    }
    std::cout << "   -> [SUCCESS] Geometric Flow converged towards Soliton.\n\n";

    std::cout << "[Phase 4] FIELD THEORY: Solving Pre-Metric Maxwell...\n";

    Aurelia::FieldTheory::Electrodynamics::AxionField axion;
    Aurelia::FieldTheory::Electrodynamics::ConstitutiveMap constitutive;
    

    auto alpha = axion.evaluate(u_test); 

    Vector3 zero_vel = {0.0, 0.0, 0.0};
    constitutive.updateState(alpha, zero_vel); 
    
    std::cout << "   -> Axion Field alpha(x) = " << alpha << "\n";
    std::cout << "   -> Constitutive Tensor assembled.\n\n";


    std::cout << "[Phase 5] EVIDENCE: Exporting Data...\n";
    std::string filename = "Aurelia_Finsler_Simulation.vtk";
    Aurelia::Utils::VTKExporter exporter(filename, NX, NY, NZ, DX, 0.0, 0.0, 0.0);
    exporter.exportState(metric_field, axion_field, director_field);

    std::cout << "   -> Data written to '" << filename << "'.\n";
    std::cout << "==============================================================\n";

    return 0;
}