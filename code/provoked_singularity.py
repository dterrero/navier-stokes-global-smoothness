# ------------------------------------------------------------------------------
# provoked_singularity.py
#
# Computes the Coherence Quotient Q(t) for a provoked singularity scenario.
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
import matplotlib.pyplot as plt
import pandas as pd

# Parameters
N = 32
L = 2 * np.pi
dt = 2e-4
nu = 0
kappa = 0
Ra = 1e5
beta = 0.0  # since kappa=0, beta will always be zero
n_steps = 3000
kc = int(N / 10)
dx = L / N

# Grid
x = np.linspace(0, L, N, endpoint=False)
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
K = fftfreq(N, d=dx) * 2 * np.pi
KX, KY, KZ = np.meshgrid(K, K, K, indexing='ij')
K2 = KX**2 + KY**2 + KZ**2
K2[K2 == 0] = 1e-10
K2_exp = K2[..., None]

# Initial fields
u = np.zeros((N, N, N, 3))
theta = 0.01 * np.sin(X) * np.sin(Y) * np.sin(Z)
u[..., 0] = 0.01 * np.sin(X) * np.cos(Y) * np.cos(Z)
u[..., 1] = -0.01 * np.cos(X) * np.sin(Y) * np.cos(Z)
u[..., 2] = 0.0

# Storage
time = np.zeros(n_steps)
Nu = np.zeros(n_steps)
Q = np.zeros(n_steps)
KE = np.zeros(n_steps)
heat_flux = np.zeros(n_steps)
mean_temp = np.zeros(n_steps)
dot_array = np.zeros(n_steps)
norm_grad_u_array = np.zeros(n_steps)
norm_A_array = np.zeros(n_steps)

# Projection
def project_div_free(u_hat):
    k_dot_u = KX * u_hat[..., 0] + KY * u_hat[..., 1] + KZ * u_hat[..., 2]
    factor = k_dot_u / K2
    u_hat[..., 0] -= KX * factor
    u_hat[..., 1] -= KY * factor
    u_hat[..., 2] -= KZ * factor
    return u_hat

# Navier-Stokes RHS
def compute_rhs(u, theta):
    grad_u = np.stack(np.gradient(u, dx, axis=(0, 1, 2)), axis=-1)
    nonlinear = np.einsum('...j,...ij->...i', u, grad_u)
    u_hat = fftn(u, axes=(0, 1, 2))
    u_hat = project_div_free(u_hat)
    viscous = ifftn(-nu * K2_exp * u_hat, axes=(0, 1, 2)).real
    buoyancy = np.zeros_like(u)
    buoyancy[..., 2] = beta * theta
    return viscous - nonlinear + buoyancy

# Q(t)
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
        Q_val = 1.0  # ✅ Perfect coherence when field is inactive
    else:
        Q_val = dot / (norm_u * norm_A + 1e-10)

    if return_diagnostics:
        return Q_val, dot, norm_u, norm_A
    else:
        return Q_val


def compute_Q_diagnostics(grad_u, kc):
    grad_A = np.zeros_like(grad_u)
    for i in range(3):
        for j in range(3):
            grad_ij_hat = fftn(grad_u[..., i, j])
            mask = (np.sqrt(KX**2 + KY**2 + KZ**2) <= kc)
            grad_A[..., i, j] = np.real(ifftn(grad_ij_hat * mask))
    dot = np.sum(grad_u * grad_A)
    norm_u = np.linalg.norm(grad_u)
    norm_A = np.linalg.norm(grad_A)
    return dot, norm_u, norm_A

# Main loop
for step in tqdm(range(n_steps)):
    if step >= 100:
        relative_step = (step - 100) / (n_steps - 100)
        force_strength = 0.1 + 2.0 * (relative_step ** 4) + 250.0 * (relative_step ** 16)

        shear = force_strength * (
            np.sin(3 * Z + 0.3 * X) * np.sin(2 * Y + 0.2 * Z) + 0.3 * np.cos(4 * Y + 0.4 * Z)
        )
        noise = np.sin(10 * X + 3 * Y) * np.cos(12 * Z + 2 * X)

        u[..., 0] += dt * force_strength * noise
        u[..., 1] += dt * 0.5 * force_strength * noise
        u[..., 2] += dt * 0.3 * force_strength * noise

        theta += force_strength * 0.01 * dt * (
            np.sin(5 * X + 0.1 * Y) * np.sin(5 * Y) * np.sin(5 * Z)
        )

    k1 = compute_rhs(u, theta)
    k2 = compute_rhs(u + 0.5 * dt * k1, theta)
    k3 = compute_rhs(u + 0.5 * dt * k2, theta)
    k4 = compute_rhs(u + dt * k3, theta)
    u += dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    # No diffusion since kappa = 0
    grad_theta = np.stack(np.gradient(theta, dx, axis=(0, 1, 2)), axis=-1)
    adv_theta = np.einsum("...i,...i->...", u, grad_theta)
    theta -= dt * adv_theta

    grad_u = np.stack(np.gradient(u, dx, axis=(0, 1, 2)), axis=-1)
    Q[step] = compute_Q(grad_u, kc)
    heat_flux[step] = np.mean(theta * u[..., 2])
    Nu[step] = 1 + heat_flux[step]  # With kappa = 0, adjust formula
    mean_temp[step] = np.mean(theta)
    time[step] = step * dt
    KE[step] = 0.5 * np.sum(u**2)
    dot_product, norm_grad_u, norm_A = compute_Q_diagnostics(grad_u, kc)
    dot_array[step] = dot_product
    norm_grad_u_array[step] = norm_grad_u
    norm_A_array[step] = norm_A

    if np.isnan(u).any() or np.isnan(theta).any():
        print(f"💥 NaN blow-up at step {step}")
        break
    if np.isnan(grad_u).any():
        print(f"💥 Blow-up in velocity gradients at step {step}")
        break

    if step % 50 == 0:
        print(f"Step {step}: Q={Q[step]:.4e}, KE={KE[step]:.4e}, Nu={Nu[step]:.2f}, "
              f"heat_flux={heat_flux[step]:.4e}, mean_temp={mean_temp[step]:.4e}")

# Save
with h5py.File("convection_gradual_blowup_Q32.h5", "w") as f:
    f.create_dataset("Q", data=Q)
    f.create_dataset("Nu", data=Nu)
    f.create_dataset("KE", data=KE)
    f.create_dataset("heat_flux", data=heat_flux)
    f.create_dataset("mean_temp", data=mean_temp)
    f.create_dataset("time_full", data=time)

# Post-analysis
alignment_angle = np.degrees(np.arccos(np.clip(Q, -1.0, 1.0)))
dot_normalized = dot_array / (norm_grad_u_array * norm_A_array + 1e-10)
dKE_dt = np.gradient(KE, time)

df = pd.DataFrame({
    "time": time,
    "Q": Q,
    "Nu": Nu,
    "KE": KE,
    "heat_flux": heat_flux,
    "mean_temp": mean_temp,
    "dot_product": dot_array,
    "norm_grad_u": norm_grad_u_array,
    "norm_A": norm_A_array,
    "alignment_angle": alignment_angle,
    "dot_normalized": dot_normalized,
    "dKE_dt": dKE_dt
})

df.to_csv("convection_gradual_blowup_Q32.csv", index=False)
print("✅ Simulation complete — diagnostics added, CSV saved.")
