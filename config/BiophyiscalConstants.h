#ifndef AURELIA_CONFIG_BIOPHYSICAL_CONSTANTS_H
#define AURELIA_CONFIG_BIOPHYSICAL_CONSTANTS_H

#include <cmath>
#include <limits>

namespace Aurelia {
namespace Config {

    // =============================================================
    // 1. MATHEMATICAL & NUMERICAL CONSTANTS
    // =============================================================
    

    constexpr double PI = 3.14159265358979323846;
    constexpr double TWO_PI = 2.0 * PI;
    constexpr double FOUR_PI = 4.0 * PI;
    constexpr double SQRT_2 = 1.41421356237309504880;

    constexpr double NUMERICAL_EPSILON = 1.0e-12; 
    constexpr double INFINITY_VAL = std::numeric_limits<double>::infinity();

    constexpr double FINITE_DIFFERENCE_STEP = 1.0e-4;
    
    // Alias for Numerical Differentiation engine (Strict Naming)
    constexpr double NUMERICAL_DIFF_STEP = FINITE_DIFFERENCE_STEP;

    /**
     * @brief Step size for Variational Derivatives (Euler-Lagrange).
     * Perturbation in function space.
     */
    constexpr double VARIATIONAL_PERTURBATION = 1.0e-5;

    /**
     * @brief Convergence tolerance for solvers.
     */
    constexpr double OPTIMIZATION_TOLERANCE = 1.0e-6;

    // =============================================================
    // 2. FUNDAMENTAL PHYSICS
    // =============================================================

    constexpr double SPEED_OF_LIGHT = 299792458.0;         // c [m/s]
    constexpr double BOLTZMANN_CONSTANT = 1.380649e-23;    // k_B [J/K]
    constexpr double H_BAR = 1.054571817e-34;              // h_bar [J.s]
    constexpr double VACUUM_PERMEABILITY = 1.25663706212e-6; // mu_0 [H/m]
    constexpr double VACUUM_PERMITTIVITY = 1.0 / (VACUUM_PERMEABILITY * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

    // =============================================================
    // 3. THERMODYNAMICS
    // =============================================================

    constexpr double TEMPERATURE_K = 293.15; // 20°C
    constexpr double THERMAL_ENERGY = BOLTZMANN_CONSTANT * TEMPERATURE_K; // k_B T
    constexpr double BETA = 1.0 / THERMAL_ENERGY;

    constexpr double ATP_ENERGY_JOULES = 8.3e-20; // ~20 k_B T

    // =============================================================
    // 4. POLYMER PHYSICS (WORM-LIKE CHAIN)
    // =============================================================

    constexpr double COLLAGEN_PERSISTENCE_LENGTH = 50.0e-9; 
    constexpr double COLLAGEN_DIAMETER = 1.5e-9;
    constexpr double COLLAGEN_CONTOUR_LENGTH = 1.0e-6;
    
    /**
     * @brief WLC Singularity Handling limits.
     * Strictly enforces valid domain for the entropic spring model.
     */
    constexpr double WLC_EXTENSION_LIMIT_RATIO = 0.99; 
    constexpr double WLC_SINGULARITY_CUTOFF    = 0.99; // 99% of Lc
    constexpr double WLC_FORCE_PENALTY         = 1.0e9;
    constexpr double WLC_ENERGY_BARRIER        = 1.0e15; // "Infinite" energy wall

    // =============================================================
    // 5. ACTIVE MATTER (MECHANO-TRANSDUCTION)
    // =============================================================

    constexpr double MESOGLEA_VISCOSITY = 1.5e-3; // Pa.s
    
    /**
     * @brief Sensitivity of cell activity to mechanical stress (Curvature).
     * Defines how strongly alpha(x) responds to R(x).
     */
    constexpr double MECHANO_SENSITIVITY = 1.0e-3; 

    /**
     * @brief Maximum scaling factor for Arrhenius activation.
     * Prevents exp() overflow in high-stress regions.
     */
    constexpr double MAX_THERMAL_SCALING = 20.0;

    /**
     * @brief Maximum biological stress limit for mesogleal cells (Pascals).
     * Used for Thermodynamic Regularization priors.
     */
    constexpr double MAX_CELL_STRESS_LIMIT = 1000.0;

    // =============================================================
    // 6. ELECTRODYNAMICS (AXION & SKEWON)
    // =============================================================

    constexpr double REFRACTIVE_INDEX_BASE = 1.338;
    constexpr double AXION_COUPLING = 1.2e-4; // g_alpha

    /**
     * @brief Non-reciprocal mixing strength (Fresnel Drag coefficient).
     * Typically very small (~v/c).
     */
    constexpr double SKEWON_MIXING_STRENGTH = 1.0e-9;

    // =============================================================
    // 7. SOLVER GRID & TIME EVOLUTION
    // =============================================================

    constexpr double GRID_DX = 1.0e-4;      // Spatial Step
    constexpr double TIME_DT = 0.01;        // Global Time Step
    constexpr int MAX_SOLVER_ITERATIONS = 5000;   

    // =============================================================
    // 8. GEODESIC INTEGRATOR SETTINGS (RK4)
    // =============================================================
    
    /**
     * @brief Parameters for the Exponential Map and Parallel Transport.
     * Ensures high-precision trajectory integration on curved manifolds.
     */
    constexpr double GEODESIC_TIME_STEP  = 0.01; // RK4 inner step
    constexpr double GEODESIC_TOTAL_TIME = 1.0;  // Normalized t=1 for Exp map

    // =============================================================
    // 9. MORPHOGENESIS & REPAIR
    // =============================================================
    
    /**
     * @brief Criterion for Symmetrization convergence.
     * When variance of curvature < Tolerance, shape is considered a Soliton.
     */
    constexpr double SYMMETRY_TOLERANCE   = 1.0e-3; 
    
    /**
     * @brief Tissue Mobility Coefficient (M).
     * Relates Chemical Potential gradient to Velocity: v = -M * grad(mu).
     */
    constexpr double TISSUE_MOBILITY      = 0.05;   

    // =============================================================
    // 10. OPTIMIZATION SOLVER (TRUST REGION)
    // =============================================================
    
    /**
     * @brief Initial radius for the Trust Region in Riemannian Newton method.
     */
    constexpr double OPTIMIZATION_INITIAL_TRUST_RADIUS = 0.1;
    
    /**
     * @brief Max inner iterations for the Steihaug-CG subproblem.
     */
    constexpr size_t      OPTIMIZATION_MAX_CG_ITER          = 20;

} 
} 

#endif 