# 2D-Advection-Solver-for-Transport-Problems
This repository contains a MATLAB implementation of a 2D transport equation solver using finite volume methods for advection. The script includes CFL-based time-stepping, velocity field visualization, and boundary condition handling, making it a robust tool for simulating transport phenomena.

# Description 
The TRANSPORT2d.m script solves 2D transport equations on a structured grid with customizable initial and boundary conditions. It uses finite volume discretization to model the advection of scalar quantities across a rectangular domain. The solver supports visualization of velocity fields, concentration distributions, and CFL-based time-stepping for stability. 

# Key features 
- CFL Condition: Ensures stable time steps based on velocity field magnitudes.
- Customizable Velocity Fields: Includes default velocity profiles with options for user-defined settings.
- Initial and Boundary Conditions: Supports Dirichlet boundary conditions and user-specified initial distributions.
- Visualization: Plots scalar concentrations and velocity field magnitudes for insight into the flow and transport processes.

# Usage
Syntax: 
```matlab
[tsteps, nsol] = TRANSPORT2d(dt0, Tend, ifshow)
```
Inputs: 
- dt0: Initial time step size. If dt0 <= 0, CFL condition determines the time step.
- Tend: End time for the simulation.
- ifshow: Visualization frequency. Set to 0 to disable visualization.

Outputs: 
- tsteps: Array of time steps.
- nsol: Numerical solution matrix containing the concentration field at the final time step.

## License
This project is licensed under the MIT License - see the LICENSE file for details.
```
Feel free to adjust any part of this README to better fit your specific needs or preferences.
