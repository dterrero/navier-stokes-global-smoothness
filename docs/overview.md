# Overview: Coherence-Based Regularity for the 3D Navier–Stokes Equations

This repository supports the proposed resolution of the [Navier–Stokes Millennium Problem](https://www.claymath.org/millennium-problems/navier%E2%80%93stokes-equation) via a novel **Coherence Decay Framework**. The solution introduces a new functional, the **Coherence Quotient** `Q(t)`, which serves as a spectral misalignment metric, enabling analytical control over regularity in 3D incompressible flow.

---

## 🔍 Key Concepts

**Coherence Quotient** `Q(t)` is defined as:

Q(t) = ||∇u(x, t) - A(x, t)||²_L²


where `A(x, t)` is a structurally aligned filtered gradient, constructed via a spectral projection onto low-frequency modes (e.g., `|k| ≤ k_c`).

### Why it matters:
- Tracks the **loss of spectral structure** in the velocity gradient.
- Offers a **new analytic route** to global smoothness — distinct from vorticity or energy-based criteria.
- Compatible with spectral simulation frameworks and verified numerically.

---

## 📁 Repository Structure

- `coherence_theory.pdf` — Full formal document with theorems and proofs
- `code/` — Python simulation scripts:
  - `simulate_NSE.py`
  - `analyze_Q_decay.py`
  - `compute_spectra.py`
- `data/` — Example datasets and visual results:
  - `results_N128.h5`
  - `plots/Q_vs_time.png`
  - `plots/vorticity.gif`
  - `plots/energy_spectrum.png`
- `docs/` — Markdown documentation (includes this `overview.md`)
- `environment.yml` — Conda environment for reproducing the results
- `CITATION.cff` — Machine-readable citation metadata

---

## 📈 Simulation Evidence

The framework is validated through simulations of the 3D Navier–Stokes equations under periodic boundary conditions using Fourier spectral methods. Highlights:

- **Resolution**: up to `N = 128³` grid points
- **Evidence**:
  - Exponential decay of `Q(t)` across all runs
  - Bounded velocity gradients `||∇u||`
  - Structural resilience under vortex stretching and stochastic perturbations

---

## 🧠 Core Mathematical Results

**Theorem 1 (Regularity from Coherence):**  
If `Q(t)` decays exponentially and the filtered field `A(x, t)` is sufficiently smooth, then the velocity field `u(x, t)` remains globally regular.

**Theorem 2 (Bounded Energy and Gradients):**  
For any smooth initial data in `H^s`, the energy remains bounded and the norm `||∇u||` remains finite under the coherence decay condition.

These results are formalized and proven in the main document: [`coherence_theory.pdf`](../coherence_theory.pdf).

---

## 🔁 Reproducibility

To run the simulations:

```bash
# Clone the repository
git clone https://github.com/dterrero/navier-stokes-global-smoothness.git
cd navier-stokes-global-smoothness

# Create environment
conda env create -f environment.yml
conda activate nse

# Run a simulation
python code/simulate_NSE.py
