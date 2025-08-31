# Overview

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![C++](https://img.shields.io/badge/language-C%2B%2B-blue)](https://isocpp.org/)
[![deal.II](https://img.shields.io/badge/deal.II-v9.4-blue)](https://dealii.org/)
[![Python](https://img.shields.io/badge/language-Python-yellow)](https://python.org/)

This project is a **proof of concept** demonstrating that **kinetic Monte Carlo (kMC)** — a method traditionally used for modeling phase change phenomena in **computational fluid dynamics (CFD)** — can be integrated into **finite element analysis (FEA)** frameworks.

The long-term goal is to benchmark this hybrid approach against classical models such as the **Johnson-Mehl-Avrami-Kolmogorov (JMAK)** model and evaluate its accuracy in simulating **microstructural evolution in thermomechanical problems**.

The motivation behind coupling kinetic Monte Carlo with FEA has been described on [my website](https://sundar.guru/projects/1_kMC_FEA/)

## Directory Structure

```
project-root/
│
├── Transient-Thermal-Base/       # ✅ Completed: FEA with temperature load mimicking additive manufacturing
│   ├── SingleElement.cpp
│   ├── temperature_input.dat
│   ├── cp_table.dat
│   ├── k_table.dat
│   ├── ...
│
├── kmc-coupled-model/              # 🔜 Upcoming: Thermal + phase change via kinetic Monte Carlo
│   └── ...
│
├── jmak-comparison-model/          # 🔜 Upcoming: Thermal + phase change via JMAK kinetics
│   └── ...
│
└── README.md
```

## Dependencies

### C++ (core simulation)

- [deal.II](https://dealii.org/) (compiled with **PETSc** and **UMFPACK** support)
- C++17 or newer
- [UMFPACK](https://people.engr.tamu.edu/davis/suitesparse.html) (part of SuiteSparse)

### Python (visualization and input generation)

- [NumPy](https://numpy.org/)
- [Matplotlib](https://matplotlib.org/)
- [Plotly](https://plotly.com/)

## Building and running

Instructions to build and run the models have been specified in the `README.md` of the respective directiory.

## License

This project is licensed under the MIT License — see the [LICENSE](./LICENSE) file for details.

## Roadmap

- [x] **Transient thermal model**
- [ ] **Thermal + phase change via kMC**
- [ ] **Thermal + phase change via JMAK**
- [ ] **Benchmark and comparison tools**

---
