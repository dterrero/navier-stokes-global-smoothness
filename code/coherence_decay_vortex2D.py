import h5py
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.ndimage import gaussian_filter

# ========== Parameters ==========
N = 128           # Grid resolution (unchanged)
T = 5000          # Number of time steps (increased from 100)
dt = 0.005        # Time step size (unchanged)
nu = 0.01         # Viscosity (unchanged)
alpha = 0.2       # Filter parameter (unchanged)
kc = 10           # Cutoff wavenumber (unchanged)

x = np.linspace(-np.pi, np.pi, N)
y = np.linspace(-np.pi, np.pi, N)
X, Y = np.meshgrid(x, y)
dx = x[1] - x[0]

# ========== Initial Conditions ==========
u = np.sin(X) * np.cos(Y) + 0.01 * np.random.randn(N, N)
v = -np.cos(X) * np.sin(Y) + 0.01 * np.random.randn(N, N)

# ========== Laplacian Function ==========
def laplacian(f, dx):
    return (
        np.roll(f, +1, axis=0) + np.roll(f, -1, axis=0) +
        np.roll(f, +1, axis=1) + np.roll(f, -1, axis=1) -
        4 * f
    ) / dx**2

# ========== HDF5 Output ==========
with h5py.File('vortex_data_128.h5', 'w') as hf:
    hf.create_dataset('x', data=x)
    hf.create_dataset('y', data=y)
    grp = hf.create_group('time_data')

    Q_int = 0
    progress = tqdm(total=T, desc='Simulating vortex dynamics')

    for step in range(T):
        # ========== Velocity Gradients ==========
        ux = np.gradient(u, x, axis=1)
        uy = np.gradient(u, y, axis=0)
        grad_sq = ux**2 + uy**2

        # ========== Filtered Structure Tensor ==========
        A = gaussian_filter(grad_sq, sigma=alpha * N / kc)

        # ========== Q and Epsilon ==========
        Q = np.mean((ux - A)**2 + (uy - A)**2)
        epsilon = nu * np.mean(ux**2 + uy**2)
        Q_int += Q * dt

        # ========== Energy ==========
        energy = np.mean(u**2 + v**2)

        if step % 10 == 0:
            print(f"Step {step:03d} | Q = {Q:.3e} | ∫Q = {Q_int:.3e} | ε = {epsilon:.3e} | "
                  f"Energy = {energy:.3e} | α·kc = {alpha * kc:.1f}, k_max = {kc * 2:.1f}")

        # ========== Save to HDF5 ==========
        step_grp = grp.create_group(f'step_{step:04d}')
        step_grp.create_dataset('u', data=u)
        step_grp.create_dataset('v', data=v)
        step_grp.attrs['Q'] = Q
        step_grp.attrs['epsilon'] = epsilon
        step_grp.attrs['energy'] = energy

        # ========== Nonlinear Terms ==========
        du_dx = np.gradient(u, x, axis=1)
        dv_dy = np.gradient(v, y, axis=0)
        u_adv = np.clip(u * du_dx, -2, 2)
        v_adv = np.clip(v * dv_dy, -2, 2)

        # ========== Update ==========
        u += dt * (nu * laplacian(u, dx) - u_adv)
        v += dt * (nu * laplacian(v, dx) - v_adv)

        progress.update(1)

    progress.close()

# ========== Vorticity Plotting Function ==========
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
    # plt.savefig(f'vortex_{"masked" if mask else "unmasked"}.png')
    plt.close()

# ========== Final Plotting ==========
with h5py.File('vortex_data_128.h5', 'r') as hf:
    last_step = hf[f'time_data/step_{T - 1:04d}']
    data = {'u': last_step['u'][...], 'v': last_step['v'][...]}

plot_vortex(data, mask=False)
plot_vortex(data, mask=True)

print("Simulation complete! Plots saved as vortex_unmasked.png and vortex_masked.png")
