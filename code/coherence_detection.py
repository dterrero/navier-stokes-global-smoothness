
import numpy as np
import h5py
from numpy.fft import fftn, ifftn, fftfreq
from tqdm import tqdm
import matplotlib.pyplot as plt

# ========== Parameters ==========
N = 64
L = 2 * np.pi
dt = 2e-4
nu = 5e-6
kappa = 5e-6
Ra = 1e4
beta = Ra * kappa / L**3
n_steps = 5000
save_every = 100
snapshot_indices = []
kc = 10  # coherence filter cutoff

# ========== Derived ==========
dx = L / N
x = np.linspace(0, L, N, endpoint=False)
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
K = fftfreq(N, d=dx) * 2 * np.pi
KX, KY, KZ = np.meshgrid(K, K, K, indexing='ij')
K2 = KX**2 + KY**2 + KZ**2
K2_nozero = np.where(K2 == 0, 1e-10, K2)
K2_exp = K2[..., None]

# ========== Initialize Fields (zero) ==========
u = np.zeros((N, N, N, 3))
theta = np.zeros((N, N, N))

# ========== Helper Functions ==========
def project_div_free(u_hat):
    k_dot_u = KX*u_hat[...,0] + KY*u_hat[...,1] + KZ*u_hat[...,2]
    factor = k_dot_u / K2_nozero
    u_hat[...,0] -= KX * factor
    u_hat[...,1] -= KY * factor
    u_hat[...,2] -= KZ * factor
    return u_hat

def check_divergence(u):
    div = (np.gradient(u[...,0], dx, axis=0) + 
           np.gradient(u[...,1], dx, axis=1) + 
           np.gradient(u[...,2], dx, axis=2))
    return np.max(np.abs(div))

def compute_rhs(u, theta):
    grad_u = np.stack(np.gradient(u, dx, axis=(0,1,2)), axis=-1)
    nonlinear = np.einsum('...j,...ij->...i', u, grad_u)
    u_hat = fftn(u, axes=(0,1,2))
    u_hat = project_div_free(u_hat)
    viscous = ifftn(-nu * K2_exp * u_hat, axes=(0,1,2)).real
    buoyancy = np.zeros_like(u)
    buoyancy[...,2] = beta * theta
    return viscous - nonlinear + buoyancy

def compute_Q(grad_u, kc):
    grad_A = np.zeros_like(grad_u)
    for i in range(3):
        for j in range(3):
            grad_ij_hat = fftn(grad_u[..., i, j])
            k_mag = np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2)
            mask = (k_mag <= kc)
            filtered_hat = grad_ij_hat * mask
            grad_A[..., i, j] = np.real(ifftn(filtered_hat))
    num = np.linalg.norm(grad_u - grad_A)
    den = np.linalg.norm(grad_u) + 1e-10
    return num / den

# ========== Diagnostics ==========
mean_temp = np.zeros(n_steps)
heat_flux = np.zeros(n_steps)
Nu = np.zeros(n_steps)
Q = np.zeros(n_steps)
time_full = []
saved_u = []
saved_theta = []

# ========== Main Loop ==========
for step in tqdm(range(n_steps)):
    t = step * dt
    time_full.append(t)

    # Start injecting both velocity and heat at step 500
    if step == 500:
        theta += 1.0 * np.sin(X) * np.sin(Y) * np.sin(Z)
        u[..., 0] = 1e-2 * np.sin(X) * np.cos(Y) * np.cos(Z)
        u[..., 1] = 1e-2 * np.cos(X) * np.sin(Y) * np.cos(Z)
        u[..., 2] = -2e-2 * np.cos(X) * np.cos(Y) * np.sin(Z)
    elif step > 500:
        theta += 0.1 * dt * np.sin(X) * np.sin(Y) * np.sin(Z)

    # Diagnostics
    mean_temp[step] = np.mean(theta)
    heat_flux[step] = np.mean(theta * u[...,2])
    Nu[step] = 1 + (L/(kappa*(np.max(theta)-np.min(theta)+1e-10)))*heat_flux[step]

    # RK4 for velocity
    k1 = compute_rhs(u, theta)
    k2 = compute_rhs(u + 0.5*dt*k1, theta)
    k3 = compute_rhs(u + 0.5*dt*k2, theta)
    k4 = compute_rhs(u + dt*k3, theta)
    u += dt/6 * (k1 + 2*k2 + 2*k3 + k4)

    # Semi-implicit temperature
    grad_theta = np.stack(np.gradient(theta, dx, axis=(0,1,2)), axis=-1)
    adv_theta = np.einsum("...i,...i->...", u, grad_theta)
    theta_hat = fftn(theta)
    theta_hat = (theta_hat*(1 - 0.5*dt*kappa*K2)) / (1 + 0.5*dt*kappa*K2)
    theta = np.real(ifftn(theta_hat)) - dt*adv_theta

    # Compute Q
    grad_u = np.stack(np.gradient(u, dx, axis=(0,1,2)), axis=-1)
    Q[step] = compute_Q(grad_u, kc) if step >= 500 else 0.0

    # Save snapshots
    if step % save_every == 0:
        saved_u.append(u.copy())
        saved_theta.append(theta.copy())
        snapshot_indices.append(step)

    # Print
    if step % 50 == 0:
        div_val = check_divergence(u)
        KE = 0.5 * np.sum(u**2)
        print(f"Step {step}: KE={KE:.4e}, Nu={Nu[step]:.2f}, ∇·u={div_val:.1e}, Q={Q[step]:.4f}")

# ========== Save Data ==========
with h5py.File("coherence_detection.h5", "w") as f:
    f.create_dataset("velocity", data=np.array(saved_u))
    f.create_dataset("temperature", data=np.array(saved_theta))
    f.create_dataset("time_full", data=np.array(time_full))
    f.create_dataset("snapshot_indices", data=np.array(snapshot_indices))
    f.create_dataset("Nu", data=Nu)
    f.create_dataset("heat_flux", data=heat_flux)
    f.create_dataset("Q", data=Q)

# ========== Plot ==========
plt.figure(figsize=(12,5))
plt.subplot(121)
plt.plot(time_full, Nu)
plt.xlabel("Time"); plt.ylabel("Nusselt Number")
plt.subplot(122)
plt.plot(time_full, heat_flux)
plt.xlabel("Time"); plt.ylabel("Heat Flux")
plt.tight_layout()
plt.savefig("diagnostics.png")
plt.close()
