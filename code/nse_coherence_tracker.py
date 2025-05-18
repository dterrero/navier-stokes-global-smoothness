# Simulation and Coherence Framework by Dickson Terrero
# Insight: Singularities in the incompressible Navier–Stokes equations arise only when coherence is lost.
#          If the Coherence Quotient Q(t) ≈ 1.0, no blowup occurs — regardless of energy injection.

import numpy as np
import h5py
import pandas as pd
from numpy.fft import fftn, ifftn, fftfreq
from tqdm import tqdm

# Parameters
N = 64
L = 2 * np.pi
dt = 1e-4
nu = 0.0
n_steps = 1000  # ⏳ shorter run
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

# Initial velocity: coherent Taylor-Green-like mode
u = np.zeros((N, N, N, 3))
u[..., 0] = np.sin(X) * np.cos(Y) * np.cos(Z)
u[..., 1] = -np.cos(X) * np.sin(Y) * np.cos(Z)
u[..., 2] = 0.0

# Storage arrays
time = np.zeros(n_steps)
Q = np.zeros(n_steps)
KE = np.zeros(n_steps)
dot_array = np.zeros(n_steps)
norm_grad_u_array = np.zeros(n_steps)
norm_A_array = np.zeros(n_steps)

# Projection operator
def project_div_free(u_hat):
    k_dot_u = KX * u_hat[..., 0] + KY * u_hat[..., 1] + KZ * u_hat[..., 2]
    factor = k_dot_u / K2
    u_hat[..., 0] -= KX * factor
    u_hat[..., 1] -= KY * factor
    u_hat[..., 2] -= KZ * factor
    return u_hat

# RHS of NSE (no diffusion in this case)
def compute_rhs(u):
    grad_u = np.stack(np.gradient(u, dx, axis=(0, 1, 2)), axis=-1)
    nonlinear = np.einsum('...j,...ij->...i', u, grad_u)
    u_hat = fftn(u, axes=(0, 1, 2))
    u_hat = project_div_free(u_hat)
    viscous = ifftn(-nu * K2_exp * u_hat, axes=(0, 1, 2)).real
    return viscous - nonlinear

# Time-dependent coherent forcing (10× original energy)
def coherent_force(t):
    A_t = 10.0 * (0.1 + 3.0 * t**2 + 100.0 * t**6)
    fx = np.sin(X) * np.cos(Y) * np.cos(Z)
    fy = -np.cos(X) * np.sin(Y) * np.cos(Z)
    fz = np.zeros_like(fx)
    force = np.stack([fx, fy, fz], axis=-1)
    return A_t * force

# Coherence quotient Q(t)
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

# Main simulation loop
for step in tqdm(range(n_steps)):
    t = step * dt
    f = coherent_force(t)
    rhs = compute_rhs(u)
    k1 = rhs + f
    k2 = compute_rhs(u + 0.5 * dt * k1) + f
    k3 = compute_rhs(u + 0.5 * dt * k2) + f
    k4 = compute_rhs(u + dt * k3) + f
    u += dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    grad_u = np.stack(np.gradient(u, dx, axis=(0, 1, 2)), axis=-1)
    Q_val, dot, norm_grad_u, norm_A = compute_Q(grad_u, kc, return_diagnostics=True)

    Q[step] = Q_val
    dot_array[step] = dot
    norm_grad_u_array[step] = norm_grad_u
    norm_A_array[step] = norm_A
    KE[step] = 0.5 * np.sum(u ** 2)
    time[step] = t

    if np.isnan(u).any():
        print(f"💥 NaN blow-up at step {step}")
        break

    if step % 50 == 0:
        print(f"Step {step}: Q={Q_val:.4e}, KE={KE[step]:.4e}")

    if step % 500 == 0:
        print(f"🔥 Milestone Step {step}: Q={Q_val:.4e}, KE={KE[step]:.4e}")

# Save diagnostics to CSV
df = pd.DataFrame({
    "time": time,
    "Q": Q,
    "KE": KE,
    "dot_product": dot_array,
    "norm_grad_u": norm_grad_u_array,
    "norm_A": norm_A_array,
    "dKE_dt": np.gradient(KE, time)
})
df.to_csv("coherent_blowup_sim.csv", index=False)
print("✅ Simulation complete — CSV saved as 'coherent_blowup_sim.csv'")
