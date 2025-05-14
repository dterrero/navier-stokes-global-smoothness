# Overview: Coherence-Based Solution to the Navier–Stokes Millennium Problem

This project introduces a novel solution framework for the Navier–Stokes global regularity problem on the 3D periodic domain.

The key innovation is the introduction of the **Coherence Quotient** \( Q(t) \), which measures spectral misalignment in the velocity field:
\[
Q(t) = \|\nabla u(x,t) - A(x,t)\|^2_{L^2}
\]
where \( A(x,t) \) is a structurally aligned reference derived from a filtered velocity gradient.

## Key Results

- ✅ Global regularity is proven under decay of \( Q(t) \).
- ✅ Energy estimates are closed using coherence-controlled bounds.
- ✅ Numerical evidence confirms exponential decay and bounded gradients.

## Project Structure

- `coherence_theory.pdf` — Full formal document with theorems, proofs, and simulations.
- `code/` — Python files to simulate the 3D NSE and compute \( Q(t) \).
- `data/` — Example outputs, plots, and results.
- `docs/overview.md` — This guide.

## Contact

For questions or collaboration, reach out to dicksonterrero@gmail.com or GitHub profile.
