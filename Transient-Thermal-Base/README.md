# Single Element Transient Thermal Model

This module implements a simple transient thermal finite element simulation for a single element using [deal.II](https://dealii.org/). It serves as the baseline for future extensions involving phase change modelling using kMC and JMAK methods.

![Temperature Profile](./temperature_profile.png)

## Features

- Time-dependent temperature boundary conditions from external input  
- Temperature-dependent material properties (thermal conductivity, specific heat)  
- Baseline for thermal benchmarking  

## Structure

- `SingleElement.cpp` – Core implementation of the thermal solver  
- `temperature_input.dat` – Input file specifying temperature vs. time  
- `cp_table.dat`, `k_table.dat` – Material property tables  
- `Solution/` – Output files (VTU, Gnuplot)  
- `generate_temperature_input.py` – Generating AM-like thermal data with empirical equations  
- `generate_pvd.py` – Generating VTU result mapping for ParaView  
- `result.py` – Plotting temperature for all nodes  
- `solution.pvd` - ParaView file for visualising results

---
