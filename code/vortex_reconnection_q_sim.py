# ------------------------------------------------------------------------------
# vortex_reconnection_q_sim.py
#
# Simulates the evolution of two anti-parallel vortex tubes in 3D space and
# computes the Coherence Quotient Q(t), Kinetic Energy (KE), and total
# gradient energy (∇u²) over time.
#
# The goal is to observe Q(t) behavior during vortex interaction and
# potential reconnection events. Coherence dynamics are recorded throughout.
#
# Author: Dickson A. Terrero
# License: CC BY-NC 4.0 — Creative Commons Attribution–NonCommercial 4.0
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

# Parameters
N = 128
L = 2 * np.pi
dx = L / N
dt = 0.01
steps = 5000
kc = 15
viscosity = 1e-3

# Grid
x = np.linspace(0, L, N, endpoint=False)
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')

# FFT wave numbers
K = fftfreq(N, d=dx) * 2 * np.pi
KX, KY, KZ = np.meshgrid(K, K, K, indexing='ij')
K2 = KX**2 + KY**2 + KZ**2
K2[0, 0, 0] = 1e-10  # Avoid division by zero

# Initial condition: Two anti-parallel vortex tubes
def vortex_tube(y0, z0, sign=1.0, radius=0.2):
    r2 = (Y - y0)**2 + (Z - z0)**2
    core = np.exp(-r2 / radius**2)
    u = np.zeros((3, N, N, N))
    u[0] = sign * core  # Only x-component
    return u

u = vortex_tube(np.pi / 2, np.pi / 2, sign=1) + vortex_tube(3 * np.pi / 2, 3 * np.pi / 2, sign=-1)

# Q(t) Coherence Function
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

    if norm_u < 1e-10 or norm_A < 1e-10:
        Q_val = 1.0  # Perfect coherence when field is inactive
    else:
        Q_val = dot / (norm_u * norm_A + 1e-10)

    if return_diagnostics:
        return Q_val, dot, norm_u, norm_A
    else:
        return Q_val

# === Save Results ===
with h5py.File("vortex_reconnection_q_results.h5", "w") as h5f:
    dset_q = h5f.create_dataset("Q", (steps,), dtype='f')
    dset_ke = h5f.create_dataset("KE", (steps,), dtype='f')
    dset_nu = h5f.create_dataset("Nu", (steps,), dtype='f')

    # Simulation Loop
    for step in tqdm(range(steps), desc="Simulating Vortex Reconnection"):
        # Compute gradient tensor
        grad_u = np.zeros((N, N, N, 3, 3))
        for i in range(3):
            for j in range(3):
                grad_u[..., i, j] = np.gradient(u[i], dx, axis=j)

        Q = compute_Q(grad_u, kc)
        KE = 0.5 * np.sum(u**2)
        Nu = np.sum(grad_u**2)

        dset_q[step] = Q
        dset_ke[step] = KE
        dset_nu[step] = Nu

        if step % 50 == 0:
            center = N // 2
            print(f"Step {step:4d} | Q(t) = {Q:.5f} | KE = {KE:.2f} | ∇u² = {Nu:.2e}")
            print("∇u @ center:\n", grad_u[center, center, center])
            print("A @ center:\n", grad_u[center, center, center])  # Optional — could replace with filtered version

        # Passive evolution (simple decay — no Navier-Stokes for now)
        for i in range(3):
            u_hat = fftn(u[i])
            u_hat *= np.exp(-viscosity * K2 * dt)
            u[i] = np.real(ifftn(u_hat))
