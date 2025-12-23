#ifndef AURELIA_CONFIG_BIOPHYSICAL_CONSTANTS_H
#define AURELIA_CONFIG_BIOPHYSICAL_CONSTANTS_H

#include <cmath>
#include <limits>

namespace Aurelia {
namespace Config {

    // =============================================================
    // 1. MATHEMATICAL & NUMERICAL CONSTANTS
    // =============================================================
    
    constexpr long double PI = 3.141592653589793238462643383279502884L;
    constexpr long double TWO_PI = 2.0L * PI;
    constexpr long double FOUR_PI = 4.0L * PI;
    constexpr long double SQRT_2 = 1.4142135623730950488016887242097L;

    constexpr long double NUMERICAL_EPSILON = 1.0e-12L; 
    constexpr long double INFINITY_VAL = std::numeric_limits<long double>::infinity();


    constexpr long double FINITE_DIFFERENCE_STEP = 1.0e-4L;
    
    // Alias for Numerical Differentiation engine (Strict Naming)
    constexpr long double NUMERICAL_DIFF_STEP = FINITE_DIFFERENCE_STEP;

    /**
     * @brief Step size for Variational Derivatives (Euler-Lagrange).
     * Perturbation in function space.
     */
    constexpr long double VARIATIONAL_PERTURBATION = 1.0e-5L;

    /**
     * @brief Convergence tolerance for solvers.
     */
    constexpr long double OPTIMIZATION_TOLERANCE = 1.0e-6L;

    // =============================================================
    // 2. FUNDAMENTAL PHYSICS
    // =============================================================

    constexpr long double SPEED_OF_LIGHT = 299792458.0L;         // c [m/s]
    constexpr long double BOLTZMANN_CONSTANT = 1.380649e-23L;    // k_B [J/K]
    constexpr long double H_BAR = 1.054571817e-34L;              // h_bar [J.s]
    constexpr long double VACUUM_PERMEABILITY = 1.25663706212e-6L; // mu_0 [H/m]
    constexpr long double VACUUM_PERMITTIVITY = 1.0L / (VACUUM_PERMEABILITY * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

    // =============================================================
    // 3. THERMODYNAMICS
    // =============================================================

    constexpr long double TEMPERATURE_K = 293.15L; // 20°C
    constexpr long double THERMAL_ENERGY = BOLTZMANN_CONSTANT * TEMPERATURE_K; // k_B T
    constexpr long double BETA = 1.0L / THERMAL_ENERGY;

    constexpr long double ATP_ENERGY_JOULES = 8.3e-20L; // ~20 k_B T

    // =============================================================
    // 4. POLYMER PHYSICS (WORM-LIKE CHAIN)
    // =============================================================

    constexpr long double COLLAGEN_PERSISTENCE_LENGTH = 50.0e-9L; 
    constexpr long double COLLAGEN_DIAMETER = 1.5e-9L;
    constexpr long double COLLAGEN_CONTOUR_LENGTH = 1.0e-6L;
    
    /**
     * @brief WLC Singularity Handling limits.
     * Strictly enforces valid domain for the entropic spring model.
     */
    constexpr long double WLC_EXTENSION_LIMIT_RATIO = 0.99L; 
    constexpr long double WLC_SINGULARITY_CUTOFF    = 0.99L; // 99% of Lc
    constexpr long double WLC_FORCE_PENALTY         = 1.0e9L;
    constexpr long double WLC_ENERGY_BARRIER        = 1.0e15L; // "Infinite" energy wall

    // =============================================================
    // 5. ACTIVE MATTER (MECHANO-TRANSDUCTION)
    // =============================================================

    constexpr long double MESOGLEA_VISCOSITY = 1.5e-3L; // Pa.s
    
    /**
     * @brief Sensitivity of cell activity to mechanical stress (Curvature).
     * Defines how strongly alpha(x) responds to R(x).
     */
    constexpr long double MECHANO_SENSITIVITY = 1.0e-3L; 

    /**
     * @brief Maximum scaling factor for Arrhenius activation.
     * Prevents exp() overflow in high-stress regions.
     */
    constexpr long double MAX_THERMAL_SCALING = 20.0L;

    /**
     * @brief Maximum biological stress limit for mesogleal cells (Pascals).
     * Used for Thermodynamic Regularization priors.
     */
    constexpr long double MAX_CELL_STRESS_LIMIT = 1000.0L;

    // =============================================================
    // 6. ELECTRODYNAMICS (AXION & SKEWON)
    // =============================================================

    constexpr long double REFRACTIVE_INDEX_BASE = 1.338L;
    constexpr long double AXION_COUPLING = 1.2e-4L; // g_alpha

    /**
     * @brief Non-reciprocal mixing strength (Fresnel Drag coefficient).
     * Typically very small (~v/c).
     */
    constexpr long double SKEWON_MIXING_STRENGTH = 1.0e-9L;

    // =============================================================
    // 7. SOLVER GRID & TIME EVOLUTION
    // =============================================================

    constexpr long double GRID_DX = 1.0e-4L;      // Spatial Step
    constexpr long double TIME_DT = 0.01L;        // Global Time Step
    constexpr int MAX_SOLVER_ITERATIONS = 5000;   

    // =============================================================
    // 8. GEODESIC INTEGRATOR SETTINGS (RK4)
    // =============================================================
    
    /**
     * @brief Parameters for the Exponential Map and Parallel Transport.
     * Ensures high-precision trajectory integration on curved manifolds.
     */
    constexpr long double GEODESIC_TIME_STEP  = 0.01L; // RK4 inner step
    constexpr long double GEODESIC_TOTAL_TIME = 1.0L;  // Normalized t=1 for Exp map

    // =============================================================
    // 9. MORPHOGENESIS & REPAIR (NEW: REMOVED MAGIC NUMBERS)
    // =============================================================
    
    /**
     * @brief Criterion for Symmetrization convergence.
     * When variance of curvature < Tolerance, shape is considered a Soliton.
     */
    constexpr long double SYMMETRY_TOLERANCE   = 1.0e-3L; 
    
    /**
     * @brief Tissue Mobility Coefficient (M).
     * Relates Chemical Potential gradient to Velocity: v = -M * grad(mu).
     */
    constexpr long double TISSUE_MOBILITY      = 0.05L;   

    // =============================================================
    // 10. OPTIMIZATION SOLVER (TRUST REGION) (NEW)
    // =============================================================
    
    /**
     * @brief Initial radius for the Trust Region in Riemannian Newton method.
     */
    constexpr long double OPTIMIZATION_INITIAL_TRUST_RADIUS = 0.1L;
    
    /**
     * @brief Max inner iterations for the Steihaug-CG subproblem.
     */
    constexpr size_t      OPTIMIZATION_MAX_CG_ITER          = 20;

} // namespace Config
} // namespace Aurelia

#endif // AURELIA_CONFIG_BIOPHYSICAL_CONSTANTS_H