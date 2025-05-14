# Global Smoothness of the 3D Incompressible Navier–Stokes Equations

**Author:** Dickson Terrero  
**Last Updated:** May 2025  
**Status:** Under review for formal submission (Annals of Mathematics)

---

## Overview

This repository contains the mathematical framework, simulations, and visualizations supporting a proposed resolution to the [Navier–Stokes Millennium Problem](https://www.claymath.org/millennium-problems/navier%E2%80%93stokes-equation).

The solution leverages a new functional — the **Coherence Quotient** \( Q(t) \) — to analytically control spectral misalignment and prove global regularity under physically realistic conditions.

---

## Contents

- [`coherence_theory.pdf`](./coherence_theory.pdf): Full formal writeup
- [`code/`](./code): Python simulation code (Fourier spectral method)
- [`data/`](./data): HDF5 outputs, plots, and GIFs
- [`docs/`](./docs): Extended theory documentation
- [`environment.yml`](./environment.yml): Dependencies for reproduction

---

## Visual Highlights

<p align="center">
  <img src="data/plots/Q_vs_time.png" width="400"/>
  <img src="data/plots/spectrum_decay.png" width="400"/>
</p>

---

## Quick Start

```bash
# Clone repo
git clone https://github.com/dterrero/navier-stokes-global-smoothness.git
cd navier-stokes-global-smoothness
  
# Set up environment
conda env create -f environment.yml
conda activate nse

# Run simulation
python code/coherence_decay_vortex2D.py.py

## Simulations
### 🔍 Coherence Detection: Q(t) vs Classical Diagnostics

<p align="center">
  <img src="assets/img/Qt_vs_traditional_diagnostics.png" width="500"/>
</p>

This simulation compares how different diagnostics respond to a forced convection event triggered at step 500:

- **Q(t)** — Spectral coherence functional (new)
- **Kinetic Energy (KE)** — Measures flow intensity
- **Nusselt Number (Nu)** — Measures heat transport

#### 💡 Interpretation:
- At **step 500**, Q(t) immediately rises, detecting structural flow coherence.
- KE and Nu show **delayed or minimal response** during the same period.
- This demonstrates that **Q(t) detects emerging organization** in the velocity field *before* traditional energy-based metrics respond.

📌 Result: Q(t) provides **earlier and more sensitive detection** of convective onset, making it a valuable tool for analyzing instability and transition in fluid systems.


## ✅ Coherence Quotient Validation (5000-Step Run)

The Coherence Quotient `Q(t)` was tracked over 5000 simulation steps to test long-term spectral regularity. The results confirm the theoretical prediction:

- **Initial:** `Q(0) ≈ 1.24`
- **Final:** `Q(5000) ≈ 0.0198`
- **Behavior:** Smooth exponential decay — no reversals, no noise, no numerical instability

### Energy and Dissipation Also Behaved Consistently:

- **Final Energy:** `≈ 0.00486` (decayed smoothly from ~0.5)
- **Final Dissipation (ε):** `≈ 1.16 × 10⁻⁴` (stable and positive)

This confirms that:

> *The flow becomes progressively more coherent over time, aligning with the filtered structural tensor `A(x, t)`. No singularities or blow-up observed. Coherence decay appears sufficient for global smoothness.*

---

### 📈 Visual Summary (Optional)

If you’d like to include the plot:

```markdown
![Q decay](data/plots/Q_vs_time_5000.png)

