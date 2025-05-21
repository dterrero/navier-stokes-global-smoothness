import numpy as np
import h5py
from numpy.fft import fftn, ifftn, fftfreq
from tqdm import tqdm

# Parameters
N = 64
L = 2 * np.pi
dx = L / N
dt = 0.01
steps = 10000
kc = 15
viscosity = 1e-3

# Grids
x = np.linspace(0, L, N, endpoint=False)
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')

# FFT wave numbers
K = fftfreq(N, d=dx) * 2 * np.pi
KX, KY, KZ = np.meshgrid(K, K, K, indexing='ij')
K2 = KX**2 + KY**2 + KZ**2
K2[0, 0, 0] = 1e-10

# Initial velocity field (shear layer)
u = np.zeros((3, N, N, N))
u[0] = np.tanh((Y - L/2) * 10)

# Q(t) computation
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
        Q_val = 1.0
    else:
        Q_val = dot / (norm_u * norm_A + 1e-10)

    if return_diagnostics:
        return Q_val, dot, norm_u, norm_A
    else:
        return Q_val

# Prepare HDF5 file
with h5py.File("shear_layer_q_results.h5", "w") as h5f:
    dset_q = h5f.create_dataset("Q", (steps,), dtype='f')
    dset_ke = h5f.create_dataset("KE", (steps,), dtype='f')
    dset_nu = h5f.create_dataset("Nu", (steps,), dtype='f')

    # Main loop
    for step in tqdm(range(steps), desc="Simulating"):
        if step == 400:
            burst = np.exp(-((X - L/2)**2 + (Y - L/2)**2 + (Z - L/2)**2) / 0.01)
            u[1] += burst * 5.0

        # Compute gradient
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
            print(f"Step {step}: Q = {Q:.4f}, KE = {KE:.2f}, Nu = {Nu:.2f}")
