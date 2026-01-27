# Aurelia Aurita

![Status](https://img.shields.io/badge/Status-Research_Prototype-blue)
![Lang](https://img.shields.io/badge/Languages-C%2B%2B17%20%7C%20VHDL-orange)

So this is the codebase for my research on Finslerian Field Theory and that Quantum Optimization stuff from the STOC '25 paper. Basically looking at *Aurelia aurita* (jellyfish) mesoglea mechanics but using some heavy math to simulate it.

I've got a dual setup here: C++ for the core logic and VHDL because I needed to offload the heavy tensor math to FPGA to get it running fast enough.

## What's actually in here

### 1. Quantum Stuff (STOC 2025)
Implemented the Unitary-Invariant Spectral Sparsification.
- **Static Envelope Sampling**: Solves the drift issue without constant updates.
- **Cubic-Regularized Newton**: Uses Lanczos iteration. Needed this for the barren plateaus problem.
- **Reverse Causality**: Tracks light-cones backwards.

### 2. Physics & Geometry
This is the Finsler geometry engine.
- **Non-Riemannian**: Computes Cartan Torsion, Chern Connection, etc. Uses the Worm-Like Chain model.
- **Field Theory**: Axion/Skewon fields from collagen structures.
- **Morphogenesis**: The actual Ricci Flow simulation for tissue regenerations.

### 3. Hardware (VHDL)
The FPGA bits.
- `RicciFlow_Streaming_Eng.vhd`: The streaming engine.
- `Chemistry`: Custom high-precision units for potentials (4th order Taylor).
- `Math`: Pipelined matrix inversion and 4th-order stencils for derivatives.

## Structure

- `src/` - The C++ stuff. Active matter logic, field theory, geometry foundation (tensors, matrices), and the theoretical solvers.
- `vhdl/` - The hardware descriptions.
- `config/` - Constants.

## How to run

You'll need GCC 9+ or Clang and OpenMP.

```bash
g++ -std=c++17 -O3 -fopenmp -I./src main.cpp -o aurelia_sim
./aurelia_sim
```

For VHDL, I use GHDL (2008 standard).

```bash
ghdl -a --std=08 vhdl/Aurelia_Types_pkg.vhd vhdl/Laplacian_Order4_Unit.vhd
ghdl -e --std=08 Stencil_Buffer_5Point
```

## Credits

**William Lee**
