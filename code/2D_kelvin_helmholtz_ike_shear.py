# ------------------------------------------------------------------------------
# 2D_kevin_helmholtz_like_shear.py
#
# Simulates a 2D Kelvin–Helmholtz–like shear flow and tracks the Coherence
# Quotient Q(t), Kinetic Energy (KE), and total gradient energy (∇u²)
# as proxies for structural alignment, energy, and dissipation.
#
# A localized coherent burst is injected at t = 400 to observe its effect
# on coherence and turbulence development in the shear layer.
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
import matplotlib.pyplot as plt
from numpy.fft import fftn, ifftn, fftfreq

# === Parameters ===
N = 192
Lx = 2 * np.pi
Ly = 2 * np.pi
dx = Lx / N
dy = Ly / N
x = np.linspace(0, Lx, N, endpoint=False)
y = np.linspace(0, Ly, N, endpoint=False)
X, Y = np.meshgrid(x, y, indexing='ij')

dt = 0.05
nu = 1e-3
n_steps = 1000
kc = 10

# === Wavenumbers ===
kx = fftfreq(N, d=dx) * 2 * np.pi
ky = fftfreq(N, d=dy) * 2 * np.pi
KX, KY = np.meshgrid(kx, ky, indexing='ij')
KZ = np.zeros_like(KX)  # 2D case

# === Initial velocity field (shear layer) ===
U0 = np.tanh((Y - Ly / 2) / 0.1)
V0 = 0.05 * np.sin(2 * np.pi * X / Lx)
u = np.stack([U0, V0, np.zeros_like(U0)], axis=-1)  # Add 3rd dim for compatibility

# === Q(t) computation function ===
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

# === Time evolution with passive dissipation ===
Q_vals = []

for step in range(n_steps):
    # Compute gradient tensor grad_u[i, j] = du_i/dx_j
    grad_u = np.zeros((N, N, 3, 3))
    for i in range(3):
        for j, spacing in enumerate([dx, dy]):
            grad_u[..., i, j] = np.gradient(u[..., i], spacing, axis=j)

    # Compute Q(t)
    Q_val = compute_Q(grad_u, kc)
    Q_vals.append(Q_val)

    # Passive decay: viscous dissipation in Fourier space
    for i in range(3):
        u_hat = fftn(u[..., i])
        u_hat *= np.exp(-nu * (KX**2 + KY**2) * dt)
        u[..., i] = np.real(ifftn(u_hat))

# === Plot Q(t) ===
plt.figure(figsize=(9, 4))
plt.plot(np.arange(n_steps) * dt, Q_vals, label='Q(t)')
plt.xlabel('Time')
plt.ylabel('Coherence Quotient Q(t)')
plt.title('Q(t) Over Time — Shear Layer Dissipation')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
