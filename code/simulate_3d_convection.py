import numpy as np
import h5py
from numpy.fft import fftn, ifftn, fftfreq
from tqdm import tqdm
import matplotlib.pyplot as plt

# ========== Parameters ==========
N = 128
L = 2 * np.pi
dt = 1e-4
nu = 5e-6
kappa = 5e-6
Ra = 1e4
beta = Ra * kappa / L**3
n_steps = 1000
save_every = 1
snapshot_indices = []

# ========== Derived Grid ==========
dx = L / N
x = np.linspace(0, L, N, endpoint=False)
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
K = fftfreq(N, d=dx) * 2 * np.pi
KX, KY, KZ = np.meshgrid(K, K, K, indexing='ij')
K2 = KX**2 + KY**2 + KZ**2
K2_nozero = np.where(K2 == 0, 1e-10, K2)
K2_exp = K2[..., None]

# ========== Adaptive Spectral Cutoff ==========
K_mag = np.sqrt(KX**2 + KY**2 + KZ**2)
k_max = np.max(K_mag)
kc = 0.6 * k_max  # Can adjust for more/less aggressive filtering

# ========== Field Initialization ==========
u = np.zeros((N, N, N, 3))
theta = np.zeros((N, N, N))

# ========== Diagnostics ==========
mean_temp = np.zeros(n_steps)
heat_flux = np.zeros(n_steps)
Nu = np.zeros(n_steps)
Q = np.zeros(n_steps)
Q_num = np.zeros(n_steps)
norm_grad_u = np.zeros(n_steps)
norm_A = np.zeros(n_steps)
KE_series = np.zeros(n_steps)

saved_u = []
saved_theta = []

# ========== Core Functions ==========
def project_div_free(u_hat):
    k_dot_u = KX * u_hat[..., 0] + KY * u_hat[..., 1] + KZ * u_hat[..., 2]
    factor = k_dot_u / K2_nozero
    u_hat[..., 0] -= KX * factor
    u_hat[..., 1] -= KY * factor
    u_hat[..., 2] -= KZ * factor
    return u_hat

def check_divergence(u):
    div = (np.gradient(u[..., 0], dx, axis=0) +
           np.gradient(u[..., 1], dx, axis=1) +
           np.gradient(u[..., 2], dx, axis=2))
    return np.max(np.abs(div))

def compute_rhs(u, theta):
    grad_u = np.stack(np.gradient(u, dx, axis=(0, 1, 2)), axis=-1)
    nonlinear = np.einsum('...j,...ij->...i', u, grad_u)
    u_hat = fftn(u, axes=(0, 1, 2))
    u_hat = project_div_free(u_hat)
    viscous = ifftn(-nu * K2_exp * u_hat, axes=(0, 1, 2)).real
    buoyancy = np.zeros_like(u)
    buoyancy[..., 2] = beta * theta
    return viscous - nonlinear + buoyancy

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

# ========== Main Simulation Loop ==========
for step in tqdm(range(n_steps)):
    # --- Inject sinusoidal + noise after step 500 ---
    if step == 499:
        theta += np.sin(X) * np.sin(Y) * np.sin(Z)
        u[..., 0] = 0.1 * np.sin(X) * np.cos(Y) * np.cos(Z)
        u[..., 1] = 0.1 * np.cos(X) * np.sin(Y) * np.cos(Z)
        u[..., 2] = -0.2 * np.cos(X) * np.cos(Y) * np.sin(Z)
    elif step >= 500:
        theta += 0.1 * dt * np.sin(X) * np.sin(Y) * np.sin(Z)
        u += 1e-3 * np.random.randn(*u.shape)  # 💥 Add subtle velocity noise

    # --- CFL check ---
    max_vel = np.max(np.linalg.norm(u, axis=-1))
    cfl = max_vel * dt / dx
    if cfl > 0.5:
        print(f"CFL warning: {cfl:.3f} at step {step}")

    # --- Diagnostics ---
    mean_temp[step] = np.mean(theta)
    heat_flux[step] = np.mean(theta * u[..., 2])
    Nu[step] = 1 + (L / (kappa * (np.max(theta) - np.min(theta) + 1e-10))) * heat_flux[step]
    KE_series[step] = 0.5 * np.mean(u ** 2)

    # --- RK4 velocity update ---
    k1 = compute_rhs(u, theta)
    k2 = compute_rhs(u + 0.5 * dt * k1, theta)
    k3 = compute_rhs(u + 0.5 * dt * k2, theta)
    k4 = compute_rhs(u + dt * k3, theta)
    u += dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    # --- Semi-implicit temperature update ---
    grad_theta = np.stack(np.gradient(theta, dx, axis=(0, 1, 2)), axis=-1)
    adv_theta = np.einsum("...i,...i->...", u, grad_theta)
    theta_hat = fftn(theta)
    theta_hat = (theta_hat * (1 - 0.5 * dt * kappa * K2)) / (1 + 0.5 * dt * kappa * K2)
    theta = np.real(ifftn(theta_hat)) - dt * adv_theta

    # --- Q(t) coherence measurement ---
    grad_u = np.stack(np.gradient(u, dx, axis=(0, 1, 2)), axis=-1)
    Q_val, dot, norm_u, norm_A_val = compute_Q(grad_u, kc, return_diagnostics=True)
    Q[step] = Q_val
    Q_num[step] = dot
    norm_grad_u[step] = norm_u
    norm_A[step] = norm_A_val

    # --- NaN check ---
    if np.isnan(u).any() or np.isnan(theta).any():
        print(f"💥 NaN detected at step {step}")
        break

    # --- Save snapshots ---
    if step % save_every == 0:
        saved_u.append(u.copy())
        saved_theta.append(theta.copy())
        snapshot_indices.append(step)

    # --- Print ---
    if step % 50 == 0:
        div_val = check_divergence(u)
        print(f"Step {step}: KE={KE_series[step]:.3e}, Nu={Nu[step]:.2f}, ∇·u={div_val:.1e}, Q={Q[step]:.4f}")


# ========== Save Data ==========
# ========== Save Data ==========
with h5py.File("convection_simulation.h5", "w") as f:
    f.create_dataset("velocity", data=np.array(saved_u))           # 1000 snapshots
    f.create_dataset("temperature", data=np.array(saved_theta))
    f.create_dataset("time", data=np.arange(n_steps) * dt)         # ✅ FULL time array
    f.create_dataset("Nu", data=Nu)
    f.create_dataset("heat_flux", data=heat_flux)
    f.create_dataset("Q", data=Q)
    f.create_dataset("Q_num", data=Q_num)
    f.create_dataset("norm_grad_u", data=norm_grad_u)
    f.create_dataset("norm_A", data=norm_A)
    f.create_dataset("KE", data=KE_series)


# ========== Plot ==========
plt.figure(figsize=(12, 5))
plt.subplot(121)
plt.plot(np.arange(n_steps) * dt, Nu)
plt.xlabel("Time"); plt.ylabel("Nusselt Number")

plt.subplot(122)
plt.plot(np.arange(n_steps) * dt, heat_flux)
plt.xlabel("Time"); plt.ylabel("Heat Flux")

plt.tight_layout()
plt.savefig("diagnostics.png")
plt.close()
