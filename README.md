
# Overview

This project is a **proof of concept** demonstrating that **kinetic Monte Carlo (kMC)** — a method traditionally used for modeling phase change phenomena in **computational fluid dynamics (CFD)** — can be integrated into **finite element analysis (FEA)** frameworks.

The long-term goal is to benchmark this hybrid approach against classical models such as the **Johnson-Mehl-Avrami-Kolmogorov (JMAK)** model and evaluate its accuracy in simulating **microstructural evolution in thermomechanical problems**.


## Directory Structure

```
project-root/
│
├── Transient-Thermal-Base/       # ✅ Completed: Single element FEA with temperature load mimicking additive manufacturing
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

* [deal.II](https://dealii.org/) (compiled with **PETSc** and **UMFPACK** support)
* C++17 or newer
* [UMFPACK](https://people.engr.tamu.edu/davis/suitesparse.html) (part of SuiteSparse)

### Python (visualization and input generation)

* [NumPy](https://numpy.org/)
* [Matplotlib](https://matplotlib.org/)
* [Plotly](https://plotly.com/)


## License

This work is licensed under the **Creative Commons Attribution-NonCommercial 4.0 International License**.
Commercial reproduction is prohibited.


## Roadmap

* [x] **Transient thermal model (single element)**
* [ ] **Thermal + phase change via kMC**
* [ ] **Thermal + phase change via JMAK**
* [ ] **Benchmark and comparison tools**

---

