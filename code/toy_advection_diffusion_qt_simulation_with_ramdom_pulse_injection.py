# ------------------------------------------------------------------------------
# toy_advection_diffusion_qt_simulation_with_randomized_pulse_injection.py
#
# Simulates a simple 3D vector field with advection–diffusion behavior,
# computes Coherence Quotient Q(t) with Fourier filtering,
# and saves results with center diagnostics.
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
#
# Please cite or link to this repository if using Q(t) in derivative works or publications:
# https://github.com/dterrero/q_collapse_mayfield
# ------------------------------------------------------------------------------

import numpy as np
import h5py
from numpy.fft import fftn, ifftn, fftfreq
from tqdm import tqdm

# Simulation parameters
N = 64
L = 2 * np.pi
dx = L / N
dt = 0.01
steps = 1000
kc = 15
viscosity = 1e-3

# Grids
x = np.linspace(0, L, N, endpoint=False)
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')

# FFT wave numbers
K = fftfreq(N, d=dx) * 2 * np.pi
KX, KY, KZ = np.meshgrid(K, K, K, indexing='ij')
K2 = KX**2 + KY**2 + KZ**2
K2[0, 0, 0] = 1e-10  # Avoid division by zero

# Initial velocity field (shear layer in x)
u = np.zeros((3, N, N, N))
u[0] = np.tanh((Y - L/2) * 10)

# Q(t) coherence function
def compute_Q(grad_u, kc, return_diagnostics=False):
    grad_A = np.zeros_like(grad_u)
    K_mag = np.sqrt(KX**2 + KY**2 + KZ**2)
    filter_mask = (K_mag <= kc)

    for i in range(3):
        for j in range(3):
            grad_ij = grad_u[..., i, j]
            grad_hat = fftn(grad_ij)
            grad_hat_filtered = grad_hat * filter_mask
            grad_A[..., i, j] = np.real(ifftn(grad_hat_filtered))

    dot = np.sum(grad_u * grad_A)
    norm_u = np.linalg.norm(grad_u)
    norm_A = np.linalg.norm(grad_A)

    Q_val = 1.0 if norm_u < 1e-10 or norm_A < 1e-10 else dot / (norm_u * norm_A + 1e-10)
    return (Q_val, dot, norm_u, norm_A) if return_diagnostics else Q_val

# Save to HDF5
with h5py.File("random_pulse_q_results.h5", "w") as h5f:
    dset_q = h5f.create_dataset("Q", (steps,), dtype='f')
    dset_ke = h5f.create_dataset("KE", (steps,), dtype='f')
    dset_nu = h5f.create_dataset("Nu", (steps,), dtype='f')

    for step in tqdm(range(steps), desc="Simulating Coherence Dynamics"):
        # Randomized pulse injection
        if step in [200, 400, 600, 800]:
            center = L/2 + (np.random.rand(3) - 0.5) * L * 0.3
            pulse = np.exp(-((X - center[0])**2 + (Y - center[1])**2 + (Z - center[2])**2) / 0.02)
            u += 2.0 * pulse * (np.random.randn(3, N, N, N))

        # Compute gradients
        grad_u = np.zeros((N, N, N, 3, 3))
        for i in range(3):
            for j in range(3):
                grad_u[..., i, j] = np.gradient(u[i], dx, axis=j)

        # Compute Q(t) and diagnostics
        Q, dot, norm_u, norm_A = compute_Q(grad_u, kc, return_diagnostics=True)
        KE = 0.5 * np.sum(u**2)
        Nu = np.sum(grad_u**2)

        dset_q[step] = Q
        dset_ke[step] = KE
        dset_nu[step] = Nu

        if step % 50 == 0:
            center_idx = N // 2
            print(f"Step {step:4d} | Q(t) = {Q:.5f} | KE = {KE:.2f} | ∇u² = {Nu:.2e}")
            print("∇u @ center:")
            print(np.round(grad_u[center_idx, center_idx, center_idx], 4))
            print("A @ center:")
            print(np.round((grad_u * (norm_A / (norm_u + 1e-10)))[center_idx, center_idx, center_idx], 4))  # Approx filtered structure




