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

<h2>Simulations</h2>

<h3>✅ Coherence Detection: Q(t) vs Classical Diagnostics</h3>

<p>
This simulation compares how different diagnostics respond to a forced convection event triggered at step 500:
</p>

<p align="center">
  <img src="assets/img/full_diagnostic_comparison_Q(s)_KE_Nu.png" width="500"/>
</p>

<ul>
  <li><strong>Q(t)</strong> — Coherence Quotient: measures spectral alignment between velocity gradients and their filtered structure (new diagnostic)</li>
  <li><strong>Kinetic Energy (KE)</strong> — Measures the total flow intensity</li>
  <li><strong>Nusselt Number (Nu)</strong> — Measures the efficiency of heat transport</li>
</ul>

<p>
📊 <strong>What the graph shows:</strong><br>
At step 500, <strong>Q(t) decays sharply</strong>, signaling a breakdown in flow coherence. In contrast, KE rises and Nu begins fluctuating only afterward.
</p>

<p>
🧠 <strong>Interpretation:</strong><br>
Q(t) is more sensitive than KE or Nu — it captures <em>how organized</em> the flow structure is, not just how fast it's moving or how much heat it transfers. 
A high KE can exist even in turbulent or incoherent states. Nu may fluctuate heavily, but it doesn't reveal structural alignment.
</p>

<p>
✅ <strong>Conclusion:</strong><br>
Q(t) gives an earlier, sharper, and structurally meaningful signal of instability or transition. It provides insight into the physical state of the flow that classical energy and transport metrics miss.
</p>

<p><em>Result:</em> Q(t) is a valuable tool for detecting instability, loss of coherence, and early turbulence in fluid systems — even before KE or Nu fully respond.</p>


<h3>✅ Coherence Quotient Validation (5000-Step Run)</h3>

<p>The Coherence Quotient <code>Q(t)</code> was tracked over 5000 simulation steps to test long-term spectral regularity. The results confirm the theoretical prediction:</p>

<ul>
  <li><strong>Initial:</strong> <code>Q(0) ≈ 1.24</code></li>
  <li><strong>Final:</strong> <code>Q(5000) ≈ 0.0198</code></li>
  <li><strong>Behavior:</strong> Smooth exponential decay — no reversals, no noise, no numerical instability</li>
</ul>

<p>Energy and Dissipation Also Behaved Consistently:</p>

<ul>
  <li><strong>Final Energy:</strong> <code>≈ 0.00486</code> (decayed smoothly from ~0.5)</li>
  <li><strong>Final Dissipation (ε):</strong> <code>≈ 1.16 × 10⁻⁴</code> (stable and positive)</li>
</ul>

<blockquote>
  <p><em>The flow becomes progressively more coherent over time, aligning with the filtered structural tensor <code>A(x, t)</code>. No singularities or blow-up observed. Coherence decay appears sufficient for global smoothness.</em></p>
</blockquote>

<hr>


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


