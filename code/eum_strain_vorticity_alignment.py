# ------------------------------------------------------------------------------
# eum_strain_vorticity_alignment.py
#
# Computes the Coherence Quotient Q(t) for 2D vortex decay and coherence diagnostics.
#
# Author: Dickson A. Terrero
# License: CC BY-NC 4.0 — Creative Commons Attribution–NonCommercial 4.0
# https://creativecommons.org/licenses/by-nc/4.0/
#
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

# ========== Parameters ==========
N = 128                    # Increased grid resolution
L = 2 * np.pi             # Domain size remains the same
dt = 2e-4                 # Reduced timestep for stability at higher N
nu = 1e-5                 # Same viscosity
n_steps = 5000           # Longer simulation
alpha = 1.0              # Coherence integral power
save_every = 100         # Save snapshots every 100 steps
snapshot_indices = [] 

# Optional noise
noise_amplitude = 1e-4
noise_start = 400
noise_end = 600

# Derived quantities
dx = L / N
K = fftfreq(N, d=dx) * 2 * np.pi
KX, KY, KZ = np.meshgrid(K, K, K, indexing='ij')
K2 = KX**2 + KY**2 + KZ**2
K2[K2 == 0] = 1e-10
K2_exp = K2[..., None]
K_mag = np.sqrt(KX**2 + KY**2 + KZ**2)

# ========== Utilities ==========

def project_div_free(u_hat):
    k_dot_u = KX * u_hat[..., 0] + KY * u_hat[..., 1] + KZ * u_hat[..., 2]
    factor = k_dot_u / K2
    u_hat[..., 0] -= factor * KX
    u_hat[..., 1] -= factor * KY
    u_hat[..., 2] -= factor * KZ
    return u_hat

def compute_rhs(u):
    grad_u = np.stack(np.gradient(u, dx, axis=(0, 1, 2)), axis=-1)
    nonlinear = np.einsum('...j,...ij->...i', u, grad_u)
    u_hat = fftn(u, axes=(0, 1, 2))
    u_hat = project_div_free(u_hat)
    viscous = ifftn(-nu * K2_exp * u_hat, axes=(0, 1, 2)).real
    return viscous - nonlinear

def init_velocity():
    u_hat = np.zeros((N, N, N, 3), dtype=np.complex128)
    decay = 1 / (1 + K2)**2
    for i in range(3):
        phase = np.exp(2j * np.pi * np.random.rand(N, N, N))
        amp = np.random.randn(N, N, N)
        u_hat[..., i] = 4.0 * decay * amp * phase
    u_hat = project_div_free(u_hat)
    return np.real(ifftn(u_hat, axes=(0, 1, 2)))

def compute_Q(grad_u, kc, return_diagnostics=False):
    grad_A = np.zeros_like(grad_u)
    filter_mask = (K_mag <= kc)
    for i in range(3):
        for j in range(3):
            grad_ij = grad_u[..., i, j]
            grad_hat = fftn(grad_ij)
            grad_hat_filtered = grad_hat * filter_mask
            grad_A[..., i, j] = np.real(ifftn(grad_hat_filtered))

    dot = np.sum(grad_u * grad_A)
    norm_u = np.linalg.norm(grad_u)
    norm_A_val = np.linalg.norm(grad_A)

    Q_val = 1.0 if norm_u < 1e-10 or norm_A_val < 1e-10 else dot / (norm_u * norm_A_val + 1e-10)

    if return_diagnostics:
        return Q_val, dot, norm_u, norm_A_val
    else:
        return Q_val

# ========== Initialization ==========
u = init_velocity()
Q = np.zeros(n_steps)
Q_int = np.zeros(n_steps)
Q_dot = np.zeros(n_steps)
Q_norm_grad_u = np.zeros(n_steps)
Q_norm_A = np.zeros(n_steps)

E = np.zeros(n_steps)
dEdt = np.zeros(n_steps)
eps = np.zeros(n_steps)
kc = np.zeros(n_steps)
vorticity_norm = np.zeros(n_steps)
Q_strain = np.zeros(n_steps)
saved_u = []

# ========== Main Loop ==========
for step in tqdm(range(n_steps)):
    if noise_start <= step < noise_end:
        u += noise_amplitude * (np.random.rand(*u.shape) - 0.5)

    grad_u = np.stack(np.gradient(u, dx, axis=(0, 1, 2)), axis=-1)

    # --- Vorticity and strain ---
    du_x = np.gradient(u[..., 0], dx, axis=(0, 1, 2))
    du_y = np.gradient(u[..., 1], dx, axis=(0, 1, 2))
    du_z = np.gradient(u[..., 2], dx, axis=(0, 1, 2))
    omega = np.stack([
        du_z[1] - du_y[2],
        du_x[2] - du_z[0],
        du_y[0] - du_x[1]
    ], axis=-1)
    omega_mag = np.sqrt(np.sum(omega**2, axis=-1))
    vorticity_norm[step] = np.mean(omega_mag)

    S = 0.5 * (grad_u + np.swapaxes(grad_u, -2, -1))
    omega_outer = np.einsum("...i,...j->...ij", omega, omega)
    A_field = np.einsum("...ij,...ij->...", S, omega_outer)
    norm_S = np.linalg.norm(S)
    norm_omega_sq = np.linalg.norm(omega)**2
    Q_strain[step] = np.mean(A_field) / (norm_S * norm_omega_sq + 1e-12)

    # --- Q(t) from spectral coherence with diagnostics ---
    eps_val = nu * np.sum(grad_u**2)
    kc_val = (eps_val / nu**3)**0.25 / (2 * np.pi) # kc_val is a diagnostic coherence scale (not used to modify the flow)
    kc_val = min(kc_val, N / 2)

    Q_val, dot, norm_u, norm_A_val = compute_Q(grad_u, kc_val, return_diagnostics=True)
    Q[step] = Q_val
    Q_dot[step] = dot
    Q_norm_grad_u[step] = norm_u
    Q_norm_A[step] = norm_A_val
    Q_int[step] = Q_val**alpha * dt + (Q_int[step - 1] if step > 0 else 0.0)

    # --- Energy and Dissipation ---
    E[step] = 0.5 * np.sum(u**2)
    if step > 0:
        dEdt[step] = -(E[step] - E[step - 1]) / dt
    eps[step] = eps_val
    kc[step] = kc_val

    if step % 100 == 0:
        print(f"Step {step}: Q={Q_val:.4e}, E={E[step]:.4e}, eps={eps_val:.4e}, Q_strain={Q_strain[step]:.4e}")

    if step % save_every == 0:
        saved_u.append(u.copy())
        snapshot_indices.append(step)

    # --- RK4 integration ---
    k1 = compute_rhs(u)
    k2 = compute_rhs(u + 0.5 * dt * k1)
    k3 = compute_rhs(u + 0.5 * dt * k2)
    k4 = compute_rhs(u + dt * k3)
    u += dt / 6 * (k1 + 2*k2 + 2*k3 + k4)

# ========== Save Output ==========
with h5py.File("eum_updated_output.h5", "w") as hf:
    hf.create_dataset("time", data=np.arange(n_steps) * dt)
    hf.create_dataset("Q", data=Q)
    hf.create_dataset("Q_integral", data=Q_int)
    hf.create_dataset("Q_dot", data=Q_dot)
    hf.create_dataset("Q_norm_grad_u", data=Q_norm_grad_u)
    hf.create_dataset("Q_norm_A", data=Q_norm_A)
    hf.create_dataset("energy", data=E)
    hf.create_dataset("dE_dt", data=dEdt)
