# Aurelia: Finslerian Active Matter Simulation

![Build Status](https://img.shields.io/badge/build-passing-brightgreen)
![License](https://img.shields.io/badge/license-MIT-blue)
![C++](https://img.shields.io/badge/std-c%2B%2B17-orange)

An HPC-optimized C++ simulation engine modeling the mesoglea of *Aurelia aurita* as a **Finslerian Manifold** driven by non-equilibrium thermodynamic flows.

## 🚀 Key Features
* **Zero-Dependency**: Written in pure C++17 with no external libraries.
* **HPC Optimized**: Uses Stack Allocation (`std::array`) and Template Metaprogramming for cache-friendly performance.
* **Geometric Rigor**: Implements Chern Connection and Cartan Torsion explicitly.
* **Paraview Export**: Native VTK output for tensor field visualization.

## 🛠️ Build Instructions

### Prerequisites
* C++17 Compliant Compiler (GCC, Clang, or MSVC)
* CMake 3.15+

### Building (Linux/macOS)
```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
./AureliaSim