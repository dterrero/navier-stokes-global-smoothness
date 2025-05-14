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
python code/simulate_NSE.py
