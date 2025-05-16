import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic

# === PARAMETERS ===
h5_path = "convection_simulation.h5"
snapshot_indices = [0, 100, 500, 999]  # Pick snapshots to analyze
L = 2 * np.pi

for snapshot_index in snapshot_indices:
    print(f"\n📊 Processing snapshot {snapshot_index}...")

    # === LOAD SNAPSHOT ===
    with h5py.File(h5_path, "r") as f:
        u = f["velocity"][snapshot_index]
        N = u.shape[0]
        dx = L / N

    # === COMPUTE WAVENUMBER GRIDS ===
    k = np.fft.fftfreq(N, d=dx) * 2 * np.pi
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    K_mag = np.sqrt(KX**2 + KY**2 + KZ**2)

    # === FOURIER TRANSFORMS ===
    u_hat = np.stack([np.fft.fftn(u[..., i]) for i in range(3)], axis=-1)

    # === ENERGY DENSITY ===
    E_density = 0.5 * np.sum(np.abs(u_hat)**2, axis=-1)

    # === VORTICITY IN FOURIER SPACE ===
    omega_hat = np.zeros_like(u_hat, dtype=complex)
    omega_hat[..., 0] = 1j * (KY * u_hat[..., 2] - KZ * u_hat[..., 1])
    omega_hat[..., 1] = 1j * (KZ * u_hat[..., 0] - KX * u_hat[..., 2])
    omega_hat[..., 2] = 1j * (KX * u_hat[..., 1] - KY * u_hat[..., 0])

    # === HELICITY DENSITY ===
    H_density = np.real(np.sum(np.conj(u_hat) * omega_hat, axis=-1))

    # === BINNING ===
    K_flat = K_mag.flatten()
    E_flat = E_density.flatten()
    H_flat = H_density.flatten()

    k_bins = np.arange(0.0, np.max(K_flat), 2 * np.pi / L)
    k_centers = 0.5 * (k_bins[:-1] + k_bins[1:])

    E_spectrum, _, _ = binned_statistic(K_flat, E_flat, bins=k_bins, statistic='mean')
    H_spectrum, _, _ = binned_statistic(K_flat, H_flat, bins=k_bins, statistic='mean')
    counts, _ = np.histogram(K_flat, bins=k_bins)

    # === NORMALIZE BY SHELL SIZE ===
    E_spectrum *= counts
    H_spectrum *= counts

    # === INTEGRAL CHECK ===
    total_energy = np.sum(E_spectrum)
    total_helicity = np.sum(np.abs(H_spectrum))
    print(f"   Total Energy ≈ {total_energy:.4e}")
    print(f"   Total Helicity ≈ {total_helicity:.4e}")

    # === SAVE TO CSV ===
    df = pd.DataFrame({
        "k": k_centers,
        "E(k)": E_spectrum,
        "H(k)": H_spectrum
    })
    csv_output_path = f"spectrum_snapshot_{snapshot_index}.csv"
    df.dropna().to_csv(csv_output_path, index=False)
    print(f"   ✅ Saved to: {csv_output_path}")

    # === PLOT ===
    plt.figure(figsize=(10, 4))

    # Energy spectrum
    plt.subplot(1, 2, 1)
    plt.loglog(k_centers, E_spectrum, label="Energy Spectrum $E(k)$")
    # Optional: add power-law reference
    inertial_range = (k_centers > 4) & (k_centers < 20)
    if np.any(inertial_range):
        ref_k = k_centers[inertial_range]
        ref_line = E_spectrum[inertial_range][0] * (ref_k / ref_k[0])**(-5/3)
        plt.loglog(ref_k, ref_line, 'k--', label=r"$k^{-5/3}$")
    plt.xlabel("Wavenumber $k$")
    plt.ylabel("$E(k)$")
    plt.grid(True)
    plt.legend()

    # Helicity spectrum
    plt.subplot(1, 2, 2)
    plt.semilogx(k_centers, H_spectrum, label="Helicity Spectrum $H(k)$", color='orange')
    plt.xlabel("Wavenumber $k$")
    plt.ylabel("$H(k)$")
    plt.grid(True)
    plt.legend()

    plt.suptitle(f"Spectrum at Snapshot {snapshot_index}", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()
