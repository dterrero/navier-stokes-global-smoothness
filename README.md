# Global Smoothness of the 3D Incompressible Navier–Stokes Equations

**Author:** Dickson Terrero & Felix A Terrero
**Last Updated:** May 2025  
**Status:**  Submitted to Annals of Mathematics and to arXiv (pending moderation), May 2025

---

<h2>Overview</h2>

<p>This repository presents a mathematical framework, simulation suite, and supporting visualizations for a proposed resolution of the 
<a href="https://www.claymath.org/millennium-problems/navier%E2%80%93stokes-equation" target="_blank">Navier–Stokes Millennium Problem</a>.</p>

<p>The approach introduces a novel diagnostic — the <strong>Coherence Quotient</strong> <code>Q(t)</code> — which quantifies spectral misalignment and serves as the foundation for a new global regularity criterion under physically realistic conditions.</p>

<blockquote>
  <p>ℹ️ <strong>How <code>Q(t)</code> is computed from the NSE:</strong><br>
  The Coherence Quotient <code>Q(t)</code> is not an external measure — it is derived directly from the Navier–Stokes velocity gradient field <code>∇u</code>.<br>
  It compares this gradient with its low-pass filtered counterpart:</p>

  <p align="center"><code>A = P<sub>k<sub>c</sub></sub> ∇u</code></p>

  <p>where <code>P<sub>k<sub>c</sub></sub></code> denotes a spectral projection onto coherent modes.</p>

  <p>The formal definition is:</p>

  <p align="center"><code>Q(t) = ⟨∇u, A⟩ / (‖∇u‖ · ‖A‖)</code></p>

  <p>This normalized inner product (cosine similarity) quantifies how well the true velocity gradient aligns with its coherent structure.<br>
  A drop in <code>Q(t)</code> indicates emerging misalignment, spectral instability, or turbulence — even before energy-based metrics detect it.</p>
</blockquote>

<hr>

<blockquote>
  <p>📝 <strong>Submission Status:</strong><br>
  This repository accompanies our formal submission to arXiv and the <em>Annals of Mathematics</em>, currently under review (May 2025).<br>
  The complete preprint — including the full theoretical framework and proof — is available in 
  <a href="./docs/Global_Smoothness_via_Coherence_Decay_in_the_3D_Navier_Stokes_Equations.pdf">Global_Smoothness_via_Coherence_Decay_in_the_3D_Navier_Stokes_Equations.pdf</a>.<br>
  While moderation is pending, this GitHub version reflects the identical content submitted for public archival and evaluation.
  </p>
</blockquote>

## Contents

- [`Global_Smoothness_via_Coherence_Decay_in_the_3D_Navier_Stokes_Equations.pdf`](./Global_Smoothness_via_Coherence_Decay_in_the_3D_Navier_Stokes_Equations.pdf) — Full formal writeup
- [`code/`](./code) — Python simulation code (Fourier spectral methods)
- [`data/`](./data) — HDF5 simulation outputs, plots, and animations
- [`docs/`](./docs) — Extended theoretical documentation
- [`environment.yml`](./environment.yml) — Reproducible conda environment

---

## Visual Highlights

<h3>✅ Coherence Detection: Q(t) vs Classical Diagnostics</h3>

<p><strong>Initial Coherence:</strong> Both simulations presented below begin with <strong>perfect spectral alignment</strong>, i.e., 
<code>Q(0) = 1.0</code> — full coherence at initialization.</p>

<p>This first simulation compares three diagnostic quantities during a forced convection event triggered at step 500:</p>

<ul>
  <li><strong>Q(t)</strong> — <em>Coherence Quotient</em>: measures spectral alignment between the full velocity gradient field and its low-pass filtered structure</li>
  <li><strong>Kinetic Energy (KE)</strong> — captures bulk flow intensity</li>
  <li><strong>Nusselt Number (Nu)</strong> — reflects convective heat transfer efficiency</li>
</ul>


<div style="display: flex; justify-content: space-between; gap: 20px;">
  <div style="flex: 1; text-align: center;">
    <img src="assets/img/plot_q_t.png" alt="Q(t) Plot" width="100%">
    <p><strong>Figure 1:</strong> Q(t) — Coherence Quotient</p>
  </div>
  <div style="flex: 1; text-align: center;">
    <img src="assets/img/plot_ke_nu_t.png" alt="KE and Nu Plot" width="100%">
    <p><strong>Figure 2:</strong> KE(t) &amp; Nu(t) — Classical Diagnostics</p>
  </div>
</div>

<p><strong>📊 What the graph shows:</strong><br>
At step 500, <code>Q(t)</code> drops sharply — signaling structural misalignment. KE and Nu respond more slowly, highlighting their limitations in capturing early instability.</p>

<p><strong>🧠 Interpretation:</strong><br>
While KE tracks energy and Nu tracks heat, only <code>Q(t)</code> reflects the internal order of the flow field. It detects breakdowns in spectral coherence well before energy-based measures do.</p>

<blockquote>
  💡 <strong>Result:</strong><br>
  <code>Q(t)</code> is a powerful structural diagnostic — capable of identifying early-stage instability, turbulence onset, and loss of smoothness far in advance of classical quantities.
</blockquote>


---

### ✅ Coherence Quotient Validation: Resolution Comparison (1000-Step Runs)

This comparison illustrates how resolution affects coherence decay over time. Two simulations — one at \( 64^3 \), the other at \( 128^3 \) — were run for 1000 steps under identical forcing.

<p align="center">
  <img src="assets/img/resolution_comparison_of_coherence_decay_128_64.png" width="500"/>
</p>

**Initial Coherence:**  
- \( Q(0) = 1.0 \) in both cases — perfect alignment at initialization

**Final Coherence at Step 1000:**  
- \( Q \approx 0.798 \) for \( 128^3 \) — faster, sharper coherence loss  
- \( Q \approx 0.850 \) for \( 64^3 \) — slower, more gradual decay

**Trend:**  
Higher resolution accelerates spectral misalignment by resolving finer-scale instabilities. The coherence framework scales naturally with grid fidelity.

**📉 Complementary Diagnostic (128³ run):**  
- Final energy: \( \approx 0.00150 \) — smoothly decayed from initial ~0.5  
- Dissipation \( \varepsilon \): *Not recorded* in this diagnostic

> 🧠 **Conclusion:**  
> \( Q(t) \) responds consistently to both resolution and physical forcing. Its decay reflects the breakdown of structural order — confirming its value as a regularity-tracking tool across scales.

---

## Quick Start

## 📄 License and Usage

This project is licensed under the **Creative Commons Attribution–NonCommercial 4.0 International License (CC BY-NC 4.0)**.

- 🔒 **Non-commercial use only** — This repository, including all code, theory, and visualizations, is intended for educational and research purposes.
- 💼 **Commercial use is not permitted** without prior written permission.
- 📚 For academic use, please cite this work using the BibTeX entry in the `NOTICE` file or via the Zenodo DOI.

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC--BY--NC%204.0-blue.svg)](https://creativecommons.org/licenses/by-nc/4.0/)


```bash
# Clone the repository
git clone https://github.com/dterrero/navier-stokes-global-smoothness.git
cd navier-stokes-global-smoothness

# Create environment
conda env create -f environment.yml
conda activate nse

# Run main 2D vortex simulation
python code/coherence_decay_vortex2D.py
