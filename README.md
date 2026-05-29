# Turbine Blade Thermal & Stress Simulator

**Live App → [turbine-thermal-stress-analysis.streamlit.app](https://turbine-thermal-stress-analysis.streamlit.app)**

A physics-based 2D steady-state heat conduction and thermal stress simulator for aeroengine turbine blade cross-sections, built in Python. Developed as an independent engineering project following manufacturing engineering work at Tata Advanced Systems Limited (Tata Centre of Excellence for Aeroengines), where real aeroengine component data is processed for clients including **GE Aerospace** and **Safran Aircraft Engines (Pratt & Whitney)**.

---

## Executive Summary

Turbine blades in high-pressure turbine (HPT) stages operate under some of the most extreme thermal conditions in engineering — hot gas temperatures exceeding 1200°C against cooling air systems maintaining blade surfaces near 400°C. Material selection and thermal management are critical to blade longevity, aerodynamic efficiency, and safety.

This simulator models the 2D steady-state temperature distribution across a NACA aerofoil cross-section using the **finite difference method (FDM)** with **Gauss-Seidel iterative relaxation** and **Robin (convective) boundary conditions**. It computes derived engineering quantities including thermal stress index, thermal diffusivity, and Biot number — and visualises heat flux vector fields across the blade interior.

---

## Features

- **7 aerospace alloys** with real thermophysical properties sourced from ASM, MatWeb, and MMPDS
- **Adjustable boundary conditions**: hot gas temperature, cooling air temperature, film cooling coefficient h
- **3 NACA aerofoil profiles**: 0012, 2412, 4412
- **Robin boundary conditions** on blade surfaces — surface temperature determined by balance of conduction (k) and convection (h), not fixed artificially
- **Heat flux vector overlay** — visualises direction and magnitude of heat flow through the blade
- **Thermal stress indexing** — σ = E·α·ΔT for all materials
- **Biot number** — quantifies convective vs conductive resistance for each material/cooling configuration
- **Material suitability warnings** — flags when T_hot approaches a material's service limit
- **Full material comparison** — side-by-side heat maps and bar charts across all 7 alloys

---

## Physics & Methodology

### Governing Equation
Steady-state 2D heat conduction (Laplace equation): ∂²T/∂x² + ∂²T/∂y² = 0
### Boundary Conditions
| Location | Condition | Type |
|---|---|---|
| Leading edge | T = T_hot (hot gas impingement) | Dirichlet |
| Trailing edge | T = T_cool | Dirichlet |
| Blade surfaces | T_surface = (k·T_interior + h·dy·T_cool) / (k + h·dy) | Robin (convective) |

### Solver
5-point finite difference stencil, Gauss-Seidel iteration, convergence tolerance 1×10⁻⁴. Typical convergence: 30–60 iterations.

### Derived Quantities
| Quantity | Formula | Physical meaning |
|---|---|---|
| Thermal diffusivity α | k / ρCp | Rate of heat propagation through material |
| Thermal stress index σ | E · α_t · ΔT | Approximate thermal fatigue stress [Pa] |
| Biot number Bi | hL / k | Convective vs conductive resistance ratio |

---

## Materials Library

| Material | Category | k (W/m·K) | T_max (°C) | Application |
|---|---|---|---|---|
| Inconel 718 | Nickel Superalloy | 11.4 | 1200 | HPT blades (GE9X, LEAP) |
| Inconel 625 | Nickel Superalloy | 9.8 | 1050 | Combustor ducts |
| Hastelloy X | Nickel Superalloy | 9.1 | 1175 | Combustion chambers |
| Ti-6Al-4V | Titanium Alloy | 7.2 | 600 | Fan blades (CFM56, GEnx) |
| Titanium CP Gr.2 | Titanium Alloy | 16.4 | 400 | Low-temp ducting |
| Stainless Steel 316L | Stainless Steel | 16.3 | 870 | LP turbine casings |
| Aluminium 6061-T6 | Aluminium Alloy | 167.0 | 300 | Nacelle structures |

---

## Project Structure
├── app.py           # Streamlit web application
├── solver.py        # FDM heat equation solver (Gauss-Seidel + Robin BCs)
├── geometry.py      # NACA 4-digit aerofoil geometry + 2D mesh generation
├── materials.py     # Aerospace alloy thermophysical property library
├── visualizer.py    # Standalone matplotlib visualisation module
└── requirements.txt

---

## Key Engineering Takeaways

- **Robin boundary conditions** produce physically accurate surface temperatures that respond to both material conductivity and cooling film intensity — unlike simplified Dirichlet (fixed-temperature) surface conditions
- **Biot number** reveals why Aluminium 6061-T6 (Bi ≈ 0.30) behaves very differently from Inconel 718 (Bi ≈ 0.044) under identical cooling: aluminium's high k means convective resistance dominates less
- **Thermal stress indexing** explains why nickel superalloys dominate HPT applications despite lower k — their lower E·α product means less thermal fatigue stress under the same ΔT

---

## Installation & Running Locally

```bash
git clone https://github.com/shaunak44oop/Turbine-Thermal-Stress-Analysis.git
cd Turbine-Thermal-Stress-Analysis
pip install -r requirements.txt
streamlit run app.py
```
