#ifndef AURELIA_UTILS_VISUALIZATION_H
#define AURELIA_UTILS_VISUALIZATION_H

#include <string>
#include <fstream>
#include <vector>
#include <iomanip>
#include <limits>

#include "../Geometry/Manifold/MetricTensor.h"
#include "../ActiveMatter/Morphogenesis/RicciFlow.h"
#include "../FieldTheory/Electrodynamics/AxionField.h"

namespace Aurelia {
namespace Utils {

    using Real = long double;
    using MetricEngine = Aurelia::Geometry::Manifold::MetricTensor;
    using PointTM = Aurelia::Geometry::Manifold::PointTM;

    class VTKExporter {
    private:
        std::string filename_;
        size_t nx_, ny_, nz_;
        Real dx_;
        Real x0_, y0_, z0_; 

    public:
    
        VTKExporter(const std::string& filename, 
                    size_t nx, size_t ny, size_t nz, 
                    Real dx,
                    Real x0 = 0.0L, Real y0 = 0.0L, Real z0 = 0.0L)
            : filename_(filename), 
              nx_(nx), ny_(ny), nz_(nz), 
              dx_(dx),
              x0_(x0), y0_(y0), z0_(z0) {}

     
        void exportState(const Aurelia::ActiveMatter::Morphogenesis::MetricFieldSoA& metric_field,
                         const std::vector<Real>& axion_field,
                         const std::vector<Real>& cartan_torsion_field,
                         const std::vector<Aurelia::Math::Vector<Real>>& director_field,
                         const std::vector<std::pair<std::string, std::vector<Real>>>& extra_scalar_fields = {}) {
            
            std::ofstream vtk(filename_);
            if (!vtk.is_open()) return;


            vtk << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 2);

 
            vtk << "# vtk DataFile Version 3.0\n";
            vtk << "Aurelia Finsler Field Theory Simulation\n";
            vtk << "ASCII\n";
            vtk << "DATASET STRUCTURED_POINTS\n";
            vtk << "DIMENSIONS " << nx_ << " " << ny_ << " " << nz_ << "\n";
            vtk << "ORIGIN " << (double)x0_ << " " << (double)y0_ << " " << (double)z0_ << "\n";
            vtk << "SPACING " << (double)dx_ << " " << (double)dx_ << " " << (double)dx_ << "\n";
            vtk << "POINT_DATA " << nx_ * ny_ * nz_ << "\n";

 
            vtk << "SCALARS Axion_Alpha double 1\n";
            vtk << "LOOKUP_TABLE default\n";
            for (const auto& val : axion_field) {
                vtk << (double)val << "\n"; 
            }

            vtk << "SCALARS Cartan_Torsion_Norm double 1\n";
            vtk << "LOOKUP_TABLE heat\n";
            for (const auto& val : cartan_torsion_field) {
                vtk << (double)val << "\n"; 
            }

            vtk << "VECTORS Collagen_Director double\n";
            for (const auto& vec : director_field) {
                vtk << (double)vec[0] << " " << (double)vec[1] << " " << (double)vec[2] << "\n";
            }

            // [SoK] Export Extra Dissonance Fields
            for (const auto& field_pair : extra_scalar_fields) {
                vtk << "SCALARS " << field_pair.first << " double 1\n";
                vtk << "LOOKUP_TABLE default\n";
                for (const auto& val : field_pair.second) {
                    vtk << (double)val << "\n";
                }
            }


            vtk << "TENSORS Finsler_Metric double\n";
            for (size_t i = 0; i < metric_field.total_size; ++i) {

                auto g = metric_field.get(i); 
                
                vtk << (double)g(0,0) << " " << (double)g(0,1) << " " << (double)g(0,2) << " "
                    << (double)g(1,0) << " " << (double)g(1,1) << " " << (double)g(1,2) << " "
                    << (double)g(2,0) << " " << (double)g(2,1) << " " << (double)g(2,2) << "\n";
            }

            vtk.close();
            std::cout << "[Viz] Exported state to " << filename_ << std::endl;
        }
    };

}
} 

#endif 