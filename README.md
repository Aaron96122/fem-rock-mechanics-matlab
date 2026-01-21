# Phase Field Fracture Simulation using Finite Element Method

This repository contains a MATLAB implementation of a phase field fracture model for simulating crack propagation in brittle materials. The code uses the finite element method (FEM) to solve the coupled system of displacement and phase field equations, employing Newton-Raphson iteration for nonlinear solution.

## Features

- **Phase Field Modeling**: Regularized crack representation using a phase field variable
- **Finite Element Analysis**: Supports 2D plane stress/strain problems with triangular and quadrilateral elements
- **Optimized Solver**: Two solver modes (optimized and unoptimized) for performance comparison
- **Visualization**: Outputs VTK files compatible with ParaView, VisIt, or similar visualization tools
- **Force-Displacement Curves**: Generates load-displacement data for post-processing

## Requirements

- **MATLAB**
## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/phase-field-fracture-fem.git
   cd phase-field-fracture-fem
   ```

2. Ensure all `.m` files and input files are in the same directory.

## Usage

### Running the Simulation

1. Open MATLAB and navigate to the repository directory.
2. Run the main program:
   ```matlab
   main_program
   ```

The simulation will:
- Read mesh data from `fract_1ca.inp`
- Perform time-stepping with Newton-Raphson iteration
- Output results to `result_1.out`, `force-disp`, and `time_*.vtk` files

### Key Parameters

Modify parameters in `main_program.m`:
- `nstep`: Total number of time steps (default: 4000)
- `nprnt`: Output interval for VTK files (default: 50)
- `isolve`: Solver mode (1 = unoptimized, 2 = optimized)
- Material properties: `constk`, `cenerg`, `constl`, `constn`, `coneta`

## Input Files

- `fract_1ca.inp`: Main input file containing mesh data, boundary conditions, and material properties
- `fract_1hca.inp`: Alternative input file (different geometry)

## Output Files

- `result_1.out`: Text file with simulation summary and computation time
- `force-disp`: Force-displacement curve data (two columns: displacement, force)
- `time_*.vtk`: VTK files for visualization at different time steps (e.g., `time_50.vtk`, `time_100.vtk`, etc.)

## Visualization

1. Install ParaView (https://www.paraview.org/) or VisIt
2. Open the VTK files to visualize:
   - Deformed mesh geometry
   - Phase field distribution (scalar field)
   - Crack propagation over time

## File Structure

```
├── main_program.m          # Main simulation driver
├── my_input_data.m         # Input file parser
├── my_initilize.m          # Initialization of variables and cracks
├── my_gauss.m              # Gauss quadrature points and weights
├── my_shape_fun.m          # Element shape functions
├── my_jacobian.m           # Jacobian and coordinate transformation
├── my_cart_deriv.m         # Pre-computed Cartesian derivatives
├── my_B_Mat.m              # Strain-displacement matrix (B)
├── bmats1.m, bmats2.m      # Vectorized B matrix functions
├── my_D_Mat.m              # Constitutive matrix (D)
├── my_stiff_v1.m           # Unoptimized stiffness assembly
├── my_stiff_v2.m           # Optimized stiffness assembly
├── my_BC_v2.m              # Boundary condition application
├── my_stress_fract_v2.m    # Stress/strain computation
├── my_residual_v2.m        # Residual force computation
├── my_vtk.m                # VTK file output
├── fract_1ca.inp           # Main input mesh file
├── fract_1hca.inp          # Alternative input mesh file
├── result_1.out            # Simulation results (generated)
├── force-disp              # Force-displacement data (generated)
└── time_*.vtk              # Visualization files (generated)
```

**Note**: The codes in this repo are the reproducible version of one project described in the book ***Porgramming Phase-Field Modeling***(https://link.springer.com/book/10.1007/978-3-319-41196-5). In this context, it is provided as-is for educational and research purposes. Please verify results for specific applications. 
