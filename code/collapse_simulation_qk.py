# ------------------------------------------------------------------------------
# collapse_simulation_qk.py
#
# Simulates a 3D vortex system with coherence-disrupting forcing and tracks:
# - Global Coherence Quotient Q(t)
# - Spectral Coherence Q(k, t) across logarithmic wavenumber bins
#
# The simulation evolves a pair of vortex tubes under mild viscosity and introduces
# random forcing during steps 100–199 to disrupt structure. The purpose is to
# observe coherence decay and how it varies across scales.
#
# No explicit structural surrender or rebirth is modeled. Q(t) and Q(k, t) are
# stored at each step for post-analysis of coherence loss and spectral collapse.
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
from tqdm import tqdm
from numpy.fft import fftn, ifftn, fftfreq

# === Parameters ===
N = 64
L = 2 * np.pi
dx = L / N
dt = 0.01
steps = 500
kc = 15
viscosity = 1e-3
n_bins = 20
output_file = "forced_collapse_qk_sim.h5"

# === Grid and FFT setup ===
x = np.linspace(0, L, N, endpoint=False)
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
kx = fftfreq(N, d=dx) * 2 * np.pi
KX, KY, KZ = np.meshgrid(kx, kx, kx, indexing='ij')
K2 = KX**2 + KY**2 + KZ**2
K2[0, 0, 0] = 1e-10
K_mag = np.sqrt(K2)

# === Initial vortex pair ===
def vortex_tube(y0, z0, sign=1.0, radius=0.2):
    r2 = (Y - y0)**2 + (Z - z0)**2
    core = np.exp(-r2 / radius**2)
    u = np.zeros((3, N, N, N))
    u[0] = sign * core
    return u

u = vortex_tube(np.pi / 2, np.pi / 2, sign=1) + vortex_tube(3 * np.pi / 2, 3 * np.pi / 2, sign=-1)

# === Utility functions ===
def compute_gradient(u):
    grad = np.zeros((N, N, N, 3, 3))
    for i in range(3):
        u_hat = fftn(u[i])
        for j, K in enumerate([KX, KY, KZ]):
            grad[..., i, j] = np.real(ifftn(1j * K * u_hat))
    return grad

def project_A(grad_u, kc):
    A = np.zeros_like(grad_u)
    mask = (K_mag <= kc)
    for i in range(3):
        for j in range(3):
            g_hat = fftn(grad_u[..., i, j])
            g_hat *= mask
            A[..., i, j] = np.real(ifftn(g_hat))
    return A

def Qk_spectrum(grad_u, A, n_bins):
    edges = np.linspace(0, np.max(K_mag), n_bins + 1)
    Qk = []
    for b in range(n_bins):
        mask = (K_mag >= edges[b]) & (K_mag < edges[b + 1])
        if not np.any(mask):
            Qk.append(np.nan)
            continue
        num = den_u = den_A = 0.0
        for i in range(3):
            for j in range(3):
                u_hat = fftn(grad_u[..., i, j])
                A_hat = fftn(A[..., i, j])
                u_shell = np.real(ifftn(u_hat * mask))
                A_shell = np.real(ifftn(A_hat * mask))
                num += np.sum(u_shell * A_shell)
                den_u += np.sum(u_shell**2)
                den_A += np.sum(A_shell**2)
        Qk.append(num / (np.sqrt(den_u * den_A) + 1e-10))
    return np.array(Qk), 0.5 * (edges[:-1] + edges[1:])

# === Forcing function ===
def coherence_disrupting_force(step):
    f = np.zeros((3, N, N, N))
    if 100 <= step < 200:  # Apply random mid-k forcing
        rng = np.random.default_rng(seed=step)
        f = 0.5 * rng.standard_normal(size=f.shape)
    return f

# === Run simulation and save results ===
with h5py.File(output_file, "w") as f:
    f.create_dataset("Q", (steps,), dtype='f4')
    f.create_dataset("Qk", (steps, n_bins), dtype='f4')
    f.create_dataset("k_centers", data=0.5 * (np.linspace(0, np.max(K_mag), n_bins + 1)[:-1] +
                                              np.linspace(0, np.max(K_mag), n_bins + 1)[1:]))

    for step in tqdm(range(steps), desc="Running Collapse Simulation"):
        grad = compute_gradient(u)
        A = project_A(grad, kc)

        dot = np.sum(grad * A)
        norm_u = np.linalg.norm(grad)
        norm_A = np.linalg.norm(A)
        Qt = dot / (norm_u * norm_A + 1e-10)
        Qk, _ = Qk_spectrum(grad, A, n_bins)

        f["Q"][step] = Qt
        f["Qk"][step] = Qk

        if step % 50 == 0:
            print(f"\nStep {step:4d} | Q(t) = {Qt:.5f}")
            print(f"Qk = {np.round(Qk, 3)}")

        # Evolve + Force
        for i in range(3):
            u_hat = fftn(u[i])
            u_hat *= np.exp(-viscosity * K2 * dt)
            u[i] = np.real(ifftn(u_hat)) + coherence_disrupting_force(step)[i]
