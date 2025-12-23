#ifndef AURELIA_ACTIVEMATTER_THERMODYNAMICS_CHEMICAL_POTENTIAL_H
#define AURELIA_ACTIVEMATTER_THERMODYNAMICS_CHEMICAL_POTENTIAL_H

#include <cmath>
#include <functional>
#include <algorithm>

#include "../../../config/BiophysicalConstants.h"
#include "../../../Geometry/Manifold/MetricTensor.h"
#include "../../../Foundation/LinearAlgebra/Matrix.h"

namespace Aurelia {
namespace ActiveMatter {
namespace Thermodynamics {

    using Real = long double;
    using PointTM = Aurelia::Geometry::Manifold::PointTM;
    using Matrix3 = Aurelia::Math::Matrix<Real>;

    class ChemicalPotential {
    private:
        Real mu_0_; 
        Real mechanosensitivity_; 
        std::function<Real(const PointTM&)> cell_density_field_;

    public:
        ChemicalPotential() {
            mu_0_ = Aurelia::Config::ATP_ENERGY_JOULES;
            mechanosensitivity_ = Aurelia::Config::MECHANO_SENSITIVITY; 
            cell_density_field_ = [](const PointTM&) { return 1.0L; };
        }

        void setCellDensity(std::function<Real(const PointTM&)> func) {
            cell_density_field_ = func;
        }

        Real compute(const PointTM& u, Real RicciScalar) const {
            Real rho = cell_density_field_(u);
            
            Real activation_energy = mechanosensitivity_ * std::abs(RicciScalar);
            Real thermal_scaling = activation_energy / Aurelia::Config::THERMAL_ENERGY;
            if (thermal_scaling > Aurelia::Config::MAX_THERMAL_SCALING) 
                thermal_scaling = Aurelia::Config::MAX_THERMAL_SCALING;

            Real activity_factor = std::exp(thermal_scaling);

            return mu_0_ * rho * activity_factor;
        }

        Matrix3 computeActiveStress(const PointTM& u, 
                                    const Matrix3& g, 
                                    Real RicciScalar) const {
            Real mu = compute(u, RicciScalar);

            return g * mu;
        }
    };

} 
} 
} 

#endif