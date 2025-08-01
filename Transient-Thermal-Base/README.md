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

## Building and running

> [!WARNING]
> This section only covers building and running with the predefined parameters. If any changes are made to the `SingleElement.cpp` file for modifying the total time or dt, the `solution.pvd` needs to be regenerated using `generate_pvd.py`

To build the files, the cleanest way is to create a build directory and run cmake from the build directory to seperate your auxillary files.

```bash
mkdir -p build && cd build
cmake ..
```

This will generate an executable file `SingleElement` in the main directory which can be run to solve the model. Following which the results can be opened on paraview.

```bash
# from inside of ./build
cd ..
./SingleElement
paraview solution.pvd
```

---
