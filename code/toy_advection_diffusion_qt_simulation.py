# ------------------------------------------------------------------------------
# toy_advection_diffusion_qt_simulation.py
#
# Simulates a simple 3D vector field with advection–diffusion behavior,
# computes Coherence Quotient Q(t) with Fourier filtering,
# and saves results with center diagnostics.
# ------------------------------------------------------------------------------

import numpy as np
import h5py
from numpy.fft import fftn, ifftn, fftfreq
from tqdm import tqdm

# Parameters
N = 64
L = 2 * np.pi
dx = L / N
dt = 0.01
steps = 1000
kc = 15
viscosity = 1e-2

# Grid setup
x = np.linspace(0, L, N, endpoint=False)
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
K = fftfreq(N, d=dx) * 2 * np.pi
KX, KY, KZ = np.meshgrid(K, K, K, indexing='ij')
K2 = KX**2 + KY**2 + KZ**2
K2[0, 0, 0] = 1e-10

# Initial velocity field (smoothed shear layer)
u = np.zeros((3, N, N, N))
u[0] = np.tanh((Y - L/2) * 10)

# Compute Q(t)
def compute_Q(grad_u, kc, return_fields=False):
    grad_A = np.zeros_like(grad_u)
    Kmag = np.sqrt(KX**2 + KY**2 + KZ**2)
    mask = (Kmag <= kc)

    for i in range(3):
        for j in range(3):
            f_hat = fftn(grad_u[..., i, j])
            f_hat_filtered = f_hat * mask
            grad_A[..., i, j] = np.real(ifftn(f_hat_filtered))

    dot = np.sum(grad_u * grad_A)
    norm_u = np.linalg.norm(grad_u)
    norm_A = np.linalg.norm(grad_A)
    Q = 1.0 if norm_u < 1e-10 or norm_A < 1e-10 else dot / (norm_u * norm_A + 1e-10)

    if return_fields:
        return Q, grad_u, grad_A
    else:
        return Q

# HDF5 setup
with h5py.File("toy_qt_advection_diffusion.h5", "w") as f:
    f.create_dataset("Q", (steps,), dtype='f')
    f.create_dataset("KE", (steps,), dtype='f')
    f.create_dataset("GradNorm", (steps,), dtype='f')

    for step in tqdm(range(steps), desc="Simulating Coherence Dynamics"):
        # Compute gradient tensor
        grad_u = np.zeros((N, N, N, 3, 3))
        for i in range(3):
            for j in range(3):
                grad_u[..., i, j] = np.gradient(u[i], dx, axis=j)

        Q, full_grad, aligned_A = compute_Q(grad_u, kc, return_fields=True)
        KE = 0.5 * np.sum(u**2)
        grad_norm = np.sum(full_grad**2)

        # Save data
        f["Q"][step] = Q
        f["KE"][step] = KE
        f["GradNorm"][step] = grad_norm

        # Print diagnostics
        if step % 50 == 0:
            print(f"Step {step:4d} | Q(t) = {Q:.5f} | KE = {KE:.2f} | ∇u² = {grad_norm:.2e}")
            center = N // 2
            print("∇u @ center:")
            print(np.round(full_grad[center, center, center], 4))
            print("A @ center:")
            print(np.round(aligned_A[center, center, center], 4))

        # Evolve u via diffusion (toy evolution)
        for i in range(3):
            u_hat = fftn(u[i])
            u_hat *= np.exp(-viscosity * K2 * dt)
            u[i] = np.real(ifftn(u_hat))
