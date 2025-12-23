
# Tugas Pengelolaan Citra S1 Teknik Informatika Universitas Bina Sarana Informatika

# Aurelia Finslerian Field Theory: A Geometric Active Matter Simulation

![Language](https://img.shields.io/badge/language-C%2B%2B17-blue.svg)
![Dependencies](https://img.shields.io/badge/dependencies-None-brightgreen.svg)
![Precision](https://img.shields.io/badge/precision-Long%20Double-orange.svg)


## 📜 Executive Summary

This repository contains a **Zero-Dependency** C++ implementation of a Unified Field Theory modeling the mesoglea of *Aurelia aurita* (Moon Jellyfish). 

Unlike traditional biomechanical models that rely on engineering approximations (Linear Elasticity, Neo-Hookean, standard Maxwell Equations), this project treats the organism as a **Finslerian Manifold** evolving under **Non-Equilibrium Thermodynamic** flows.

**Key Innovations:**
1.  **Microscopic Derivation:** The metric tensor $g_{ij}(x,y)$ is derived from the statistical mechanics of **Worm-Like Chains (WLC)**, not phenomenological curve-fitting.
2.  **Rigorous Geometry:** Implements the **Chern Connection** and calculates **Cartan Torsion** to quantify non-Riemannian nonlinearity.
3.  **Chiral Electrodynamics:** Replaces standard Maxwell theory with **Pre-Metric Electrodynamics**, incorporating **Axion Fields** ($\alpha$) to model the topological magnetoelectric effect of the chiral collagen network.
4.  **Active Morphogenesis:** Models tissue repair and symmetrization using **Diffusive-Reactive Ricci Flow**, driven by cellular chemical potentials ($\mu_{bio}$).

---

## 🏛️ Addressing Theoretical Critiques

This codebase was architected specifically to resolve fundamental flaws identified in previous "Latex-like" models:

### 1. Mathematical Rigor (The Geometry)
* **Critique:** "Using a State Vector simplifies the Slit Tangent Bundle."
* **Resolution:** We implement the full **Slit Tangent Bundle** structure $\mathring{T}M = TM \setminus \{0\}$. The metric is dynamic and direction-dependent. We compute the **Cartan Torsion Tensor** $C_{ijk}$ explicitly. If $||C_{ijk}|| \neq 0$, the space is proven to be Finslerian.
* **Core Module:** `src/Geometry/Manifold/`

### 2. Physical Foundations (The Fields)
* **Critique:** "Standard Maxwell equations ignore the chirality of collagen."
* **Resolution:** We implement **Hehl-Obukhov Pre-Metric Electrodynamics**. The Constitutive Tensor $\chi^{abcd}$ (36 components) includes the **Axion Field** (Pseudoscalar) and **Skewon Field** (Non-reciprocal mixing), derived from the helical symmetry of Type-0 collagen.
* **Core Module:** `src/FieldTheory/Electrodynamics/`

### 3. Biological Agency (The Thermodynamics)
* **Critique:** "The model treats tissue as dead rubber."
* **Resolution:** We treat the tissue as **Active Matter**. The metric evolves via **Ricci Flow** coupled to ATP hydrolysis energy ($\mu_{bio}$). The **Entropy Production Rate** is monitored to ensure thermodynamic consistency (Second Law).
* **Core Module:** `src/ActiveMatter/`

---

## 📂 Project Structure

The architecture reflects the ontological hierarchy of the theory:

```text
Aurelia_Finsler_FieldTheory/
├── config/                  # Universal Biophysical Constants (CODATA/WLC)
├── src/
│   ├── Foundation/          # MANUAL MATH ENGINE (No External Libraries)
│   │   ├── LinearAlgebra/   # Matrix, Tensor4, Complex, Eigensolver
│   │   └── Calculus/        # NumericalDiff (4th Order), Gauss-Legendre Integration
│   │
│   ├── StatisticalMech/     # FROM POLYMER TO METRIC
│   │   ├── Polymer/         # Worm-Like Chain (WLC) Physics
│   │   └── Microstructure/  # Orientation Distribution Function (ODF)
│   │
│   ├── Geometry/            # DIFFERENTIAL GEOMETRY
│   │   ├── Manifold/        # Metric Tensor g_ij, Cartan Torsion C_ijk
│   │   ├── Connection/      # Chern Connection Gamma^i_jk
│   │   └── Invariants/      # Flag Curvature, Geodesic Flow (Exp Map)
│   │
│   ├── FieldTheory/         # ELECTRODYNAMICS & LAGRANGIAN
│   │   ├── Electrodynamics/ # Axion Field, Pre-Metric Maxwell
│   │   └── Lagrangian/      # Action Functional S, Euler-Lagrange Solver
│   │
│   ├── ActiveMatter/        # BIOLOGY & MORPHOGENESIS
│   │   ├── Thermodynamics/  # Chemical Potential, Dissipation
│   │   └── Morphogenesis/   # Ricci Flow, Symmetrization
│   │
│   ├── Inversion/           # OPTIMIZATION ALGORITHMS
│   │   └── RiemannianOpt/   # Manifold Newton-CG, Parallel Transport
│   │
│   └── Utils/               # VISUALIZATION & IO
│       └── Visualization.h  # VTK Exporter for Paraview
│
└── main.cpp                 # Simulation Orchestrator