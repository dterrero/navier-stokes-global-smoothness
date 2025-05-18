# ------------------------------------------------------------------------------
# coherence_decay_vortex2D.py
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



import h5py
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.ndimage import gaussian_filter

# ========== Parameters ==========
N = 128
T = 5000
dt = 0.005
nu = 0.01
alpha = 0.2
kc = 10

x = np.linspace(-np.pi, np.pi, N)
y = np.linspace(-np.pi, np.pi, N)
X, Y = np.meshgrid(x, y)
dx = x[1] - x[0]

# ========== Initial Conditions ==========
u = np.sin(X) * np.cos(Y) + 0.01 * np.random.randn(N, N)
v = -np.cos(X) * np.sin(Y) + 0.01 * np.random.randn(N, N)

# ========== Laplacian ==========
def laplacian(f, dx):
    return (
        np.roll(f, +1, axis=0) + np.roll(f, -1, axis=0) +
        np.roll(f, +1, axis=1) + np.roll(f, -1, axis=1) -
        4 * f
    ) / dx**2

# ========== HDF5 Save ==========
with h5py.File('vortex_data_128.h5', 'w') as hf:
    hf.create_dataset('x', data=x)
    hf.create_dataset('y', data=y)
    grp = hf.create_group('time_data')

    Q_int = 0
    progress = tqdm(total=T, desc='Simulating vortex dynamics')

    for step in range(T):
        # === Compute Gradients ===
        ux = np.gradient(u, dx, axis=1)
        uy = np.gradient(u, dx, axis=0)

        # === Filter Each Component ===
        Ax = gaussian_filter(ux, sigma=alpha * N / kc)
        Ay = gaussian_filter(uy, sigma=alpha * N / kc)

        # === Coherence Quotient Q(t) ===
        dot = ux * Ax + uy * Ay
        norm_grad = np.sqrt(ux**2 + uy**2)
        norm_A = np.sqrt(Ax**2 + Ay**2)

        Q = np.sum(dot) / (np.linalg.norm(norm_grad) * np.linalg.norm(norm_A) + 1e-10)

        # === Dissipation Rate ===
        epsilon = nu * np.mean(ux**2 + uy**2)

        # === Accumulate ∫Q dt ===
        Q_int += Q * dt

        # === Energy ===
        energy = np.mean(u**2 + v**2)

        # === Logging ===
        if step % 10 == 0:
            print(f"Step {step:04d} | Q = {Q:.4f} | ∫Q = {Q_int:.4e} | ε = {epsilon:.4e} | E = {energy:.4f}")

        # === Save ===
        step_grp = grp.create_group(f'step_{step:04d}')
        step_grp.create_dataset('u', data=u)
        step_grp.create_dataset('v', data=v)
        step_grp.attrs['Q'] = Q
        step_grp.attrs['epsilon'] = epsilon
        step_grp.attrs['energy'] = energy

        # === Advective Terms ===
        du_dx = np.gradient(u, dx, axis=1)
        dv_dy = np.gradient(v, dx, axis=0)
        u_adv = np.clip(u * du_dx, -2, 2)
        v_adv = np.clip(v * dv_dy, -2, 2)

        # === Update ===
        u += dt * (nu * laplacian(u, dx) - u_adv)
        v += dt * (nu * laplacian(v, dx) - v_adv)

        progress.update(1)

    progress.close()

# ========== Vorticity Plotting ==========
def plot_vortex(data, mask=False):
    u = data['u']
    v = data['v']
    vort = np.gradient(v, axis=1) - np.gradient(u, axis=0)

    if mask:
        kx = np.fft.fftfreq(N, x[1] - x[0]) * 2 * np.pi
        ky = np.fft.fftfreq(N, y[1] - y[0]) * 2 * np.pi
        KX, KY = np.meshgrid(kx, ky, indexing='ij')
        mask_array = (np.sqrt(KX ** 2 + KY ** 2) < kc).astype(float)
        vort = np.fft.ifft2(np.fft.fft2(vort) * mask_array).real

    plt.figure(figsize=(10, 8))
    plt.imshow(vort, cmap='viridis', extent=[-np.pi, np.pi, -np.pi, np.pi])
    plt.colorbar(label='Vorticity')
    plt.quiver(X[::8, ::8], Y[::8, ::8], u[::8, ::8], v[::8, ::8], color='white')
    plt.title(f'Vortex Structure {"(Filtered)" if mask else "(Unfiltered)"}')
    plt.savefig(f'vortex_{"masked" if mask else "unmasked"}.png')
    plt.close()

# ========== Final Plot ==========
with h5py.File('vortex_data_128.h5', 'r') as hf:
    last_step = hf[f'time_data/step_{T - 1:04d}']
    data = {'u': last_step['u'][...], 'v': last_step['v'][...]}

plot_vortex(data, mask=False)
plot_vortex(data, mask=True)

print("✅ Simulation complete! Final plots saved.")
