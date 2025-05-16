import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fftfreq, fftn, ifftn

# --- Parameters ---
N = 128           # Grid size
L = 2 * np.pi     # Domain length
dx = L / N

# --- Auto-compute spectral grid ---
K = fftfreq(N, d=dx) * 2 * np.pi
KX, KY, KZ = np.meshgrid(K, K, K, indexing='ij')
K_mag = np.sqrt(KX**2 + KY**2 + KZ**2)
k_max = np.max(np.abs(K))  # Should be N/2 * 2π/L = 64

# --- Recommended cutoffs ---
kc = 0.6 * k_max                # Filter cutoff (~60% of max)
kc_dealias = (2 / 3) * k_max   # Standard 2/3 rule

# --- Print diagnostics ---
print(f"Spectral grid max k: {k_max:.2f}")
print(f"Filter cutoff kc: {kc:.2f}")
print(f"Dealiasing cutoff kc_dealias: {kc_dealias:.2f}")

# --- Generate smooth random test field ---
np.random.seed(0)
u = np.random.randn(3, N, N, N) * np.exp(-K_mag**2 / (2 * kc**2))

# --- Compute ∇u and transform ---
grad_u = np.gradient(u, dx, axis=(1, 2, 3))
grad_u_stack = np.stack(grad_u)  # shape: (3, 3, N, N, N)
grad_u_hat = fftn(grad_u_stack, axes=(2, 3, 4))

# --- Apply filter mask (|K| ≤ kc) ---
filter_mask = (K_mag <= kc)
filter_mask_4d = filter_mask[None, None, :, :, :]  # Expand to (3, 3, N, N, N)
Pkc_grad_u_hat = grad_u_hat * filter_mask_4d
A = ifftn(Pkc_grad_u_hat, axes=(2, 3, 4)).real

# --- Filter energy diagnostic ---
removed_energy = np.linalg.norm((1 - filter_mask_4d) * grad_u_hat)**2
total_energy = np.linalg.norm(grad_u_hat)**2
print(f"[Filter Check] High-mode energy removed: {removed_energy / total_energy:.2%}")

# --- Dealiasing mask (2/3 rule) ---
dealias_mask = (
    (np.abs(KX) < kc_dealias) &
    (np.abs(KY) < kc_dealias) &
    (np.abs(KZ) < kc_dealias)
)
print(f"[Dealiasing Check] Modes preserved: {np.sum(dealias_mask)} / {N**3} = {100 * np.sum(dealias_mask) / N**3:.2f}%")

# --- Visualization ---
mid = N // 2
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.title("Filter Mask (|K| ≤ kc)")
plt.imshow(K_mag[:, :, mid] <= kc, cmap='viridis', origin='lower')
plt.colorbar()

plt.subplot(1, 2, 2)
plt.title("Dealiasing Mask (2/3 Rule)")
plt.imshow(dealias_mask[:, :, mid], cmap='viridis', origin='lower')
plt.colorbar()

plt.tight_layout()
plt.show()
