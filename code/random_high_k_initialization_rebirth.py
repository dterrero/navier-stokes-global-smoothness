# ------------------------------------------------------------------------------
# random_high_k_initialization_rebirth.py
#
# Simulates the evolution of coherence and energy in a 3D fluid flow
# initialized with random high-wavenumber (high-k) energy content.
# Computes the Coherence Quotient Q(t), total energy E(t), and
# gradient supremum norm ||∇u||∞ (as a proxy for dissipation intensity)
# over 1500 steps. The simulation tracks whether coherent structure
# can emerge (“rebirth”) from an initially disordered, small-scale-dominated field.

#
# Author: Dickson A. Terrero
# License: CC BY-NC 4.0
# https://creativecommons.org/licenses/by-nc/4.0/
# ------------------------------------------------------------------------------
# ⚠️ Usage Notice:
# This script and associated formulas are shared for **educational and research purposes only**.
# Commercial use is **not permitted** under the terms of the license.
#
# If you wish to use this method in a commercial product or service,
# please contact the author to discuss licensing terms.
# ------------------------------------------------------------------------------

import numpy as np
import h5py
from tqdm import tqdm
from numpy.fft import fftn, ifftn, fftfreq

# === Parameters ===
N = 128
L = 2 * np.pi
dx = L / N
dt = 0.01
steps = 1500
kc = 15
viscosity = 1e-3
output_file = "simulation_2_random_highk.h5"

# === Grid and FFT setup ===
x = np.linspace(0, L, N, endpoint=False)
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
kx = fftfreq(N, d=dx) * 2 * np.pi
KX, KY, KZ = np.meshgrid(kx, kx, kx, indexing='ij')
K2 = KX**2 + KY**2 + KZ**2
K2[0, 0, 0] = 1e-10
K_mag = np.sqrt(K2)

# === Initial condition: high-k random noise ===
rng = np.random.default_rng(seed=42)
u = rng.standard_normal((3, N, N, N)) * 0.05
for i in range(3):
    u_hat = fftn(u[i])
    u_hat *= (K_mag > kc)  # Keep only high-k
    u[i] = np.real(ifftn(u_hat))

# === Gradient computation ===
def compute_gradient(u):
    grad = np.zeros((N, N, N, 3, 3))
    for i in range(3):
        u_hat = fftn(u[i])
        for j, K in enumerate([KX, KY, KZ]):
            grad[..., i, j] = np.real(ifftn(1j * K * u_hat))
    return grad

# === Coherence projection & Q(t) ===
def compute_Q(grad_u, kc, return_A=False):
    grad_A = np.zeros_like(grad_u)
    filter_mask = (KX**2 + KY**2 + KZ**2 <= kc**2)

    for i in range(3):
        for j in range(3):
            grad_ij = grad_u[..., i, j]
            grad_hat = fftn(grad_ij)
            grad_hat_filtered = grad_hat * filter_mask
            grad_A[..., i, j] = np.real(ifftn(grad_hat_filtered))

    dot = np.sum(grad_u * grad_A)
    norm_u = np.linalg.norm(grad_u)
    norm_A = np.linalg.norm(grad_A)

    if norm_u < 1e-10:
        Q_val = 1.0  # No motion, trivially coherent
    elif norm_A < 1e-10:
        Q_val = 0.0  # No coherent structure to align with
    else:
        Q_val = dot / (norm_u * norm_A + 1e-10)

    return (Q_val, grad_A) if return_A else Q_val


# === Forcing function: small noise between steps 100–200 ===
def coherence_disrupting_force(step):
    f = np.zeros((3, N, N, N))
    if 100 <= step < 200:
        rng = np.random.default_rng(seed=step)
        f = 0.05 * rng.standard_normal(size=f.shape)
    return f

# === Run Simulation ===
with h5py.File(output_file, "w") as f:
    f.create_dataset("Q", (steps,), dtype='f4')
    f.create_dataset("E", (steps,), dtype='f4')
    f.create_dataset("grad_Linf", (steps,), dtype='f4')

    for step in tqdm(range(steps), desc="Running [random_highk]"):
        grad = compute_gradient(u)

        if step < 300:
            Qval, A = compute_Q(grad, kc, return_A=True)
        elif step < 400:
            Qval, A = 0.0, np.zeros_like(grad)
        else:
            Qval, A = compute_Q(grad, kc, return_A=True)

        E = 0.5 * np.sum(u**2)
        grad_max = np.max(np.abs(grad))

        f["Q"][step] = Qval
        f["E"][step] = E
        f["grad_Linf"][step] = grad_max

        if step % 50 == 0:
            print(f"\nStep {step:4d} | Q(t) = {Qval:.5f} | E = {E:.4f} | ∇u_L∞ = {grad_max:.3f}")
            print("Sample ∇u[0,0,0] =", np.round(grad[0, 0, 0], 3))
            print("Sample A[0,0,0]  =", np.round(A[0, 0, 0], 3))
            if Qval < 0.05:
                print("⚠️  Coherence collapse: Q(t) < 0.05")
            if grad_max > 100:
                print("⚠️  Extreme gradients: ∇u_L∞ > 100")

        # Evolve velocity field
        f_force = coherence_disrupting_force(step)
        for i in range(3):
            u_hat = fftn(u[i])
            u_hat *= np.exp(-viscosity * K2 * dt)
            u[i] = np.real(ifftn(u_hat)) + f_force[i]
