# ------------------------------------------------------------------------------
# run_q_time_series_with_artifacts_removal.py
#
# Computes the Coherence Quotient Q(t) using velocity gradient data from
# HRRR reanalysis (Herbie Main). This implementation reflects the theoretical
# structure described in EUM-based turbulence diagnostics.
#
# Author: Dickson A. Terrero
# License: CC BY-NC 4.0 — Creative Commons Attribution–NonCommercial 4.0
# https://creativecommons.org/licenses/by-nc/4.0/
#
# ----------------------------------------------------------------------
# ⚠️ Usage Notice:
# This script and associated formulas are shared for **educational and
# research purposes only**. Commercial use is **not permitted** under the
# terms of the license.
#
# If you wish to use this method in a commercial setting or product,
# please contact the author to discuss licensing terms.
#
# Cite or link to this repository if using Q(t) in derivative works or publications:
# https://github.com/dterrero/navier-stokes-global-smoothness/tree/main/herbie_qt_collapse
# ------------------------------------------------------------------------------

import numpy as np
from numpy.fft import fft2, ifft2, fftfreq
from scipy.interpolate import RegularGridInterpolator

def run_q_diagnostic_cleaned(ds, time_exact, level, lat_min, lat_max, lon_min, lon_max, kc):
    try:
        ds_slice = ds.sel(time=np.datetime64(time_exact))
    except KeyError:
        print(f"⚠️ Time {time_exact} not found in dataset.")
        return

    # Choose level with tighter control
    ds_slice = ds_slice.sel(isobaricInhPa=level, method="nearest")

    u = ds_slice['u'].values
    v = ds_slice['v'].values
    lat = ds_slice['latitude'].values
    lon = ds_slice['longitude'].values

    npts = 64
    lat_new = np.linspace(lat_min, lat_max, npts + 8)[4:-4]  # add buffer to reduce edge effects
    lon_new = np.linspace(lon_min, lon_max, npts + 8)[4:-4]
    lat_grid, lon_grid = np.meshgrid(lat_new, lon_new, indexing='ij')
    points = np.column_stack((lat_grid.ravel(), lon_grid.ravel()))

    u_interp = RegularGridInterpolator((lat[:, 0], lon[0, :]), u, bounds_error=False, fill_value=np.nan)
    v_interp = RegularGridInterpolator((lat[:, 0], lon[0, :]), v, bounds_error=False, fill_value=np.nan)

    u_grid = u_interp(points).reshape((npts, npts))
    v_grid = v_interp(points).reshape((npts, npts))

    # Mask NaNs introduced by boundary gaps
    mask_valid = ~np.isnan(u_grid) & ~np.isnan(v_grid)
    if not np.all(mask_valid):
        print("⚠️ Skipping step due to missing values after interpolation.")
        return

    vel = np.stack([u_grid, v_grid], axis=-1)
    dy = np.nanmean(np.abs(np.diff(lat_grid, axis=0)))
    dx = np.nanmean(np.abs(np.diff(lon_grid, axis=1)))

    grad_u = np.zeros((npts, npts, 2, 2))
    for i in range(2):
        for j, spacing in enumerate([dy, dx]):
            grad_u[..., i, j] = np.gradient(vel[..., i], spacing, axis=j)

    # === Filtering high-frequency filamentary noise ===
    KX = fftfreq(npts, d=dx) * 2 * np.pi
    KY = fftfreq(npts, d=dy) * 2 * np.pi
    KX, KY = np.meshgrid(KX, KY)
    Kmag = np.sqrt(KX**2 + KY**2)
    filter_mask = (Kmag <= kc)

    grad_A = np.zeros_like(grad_u)
    for i in range(2):
        for j in range(2):
            f_hat = fft2(grad_u[..., i, j])
            f_hat_filtered = f_hat * filter_mask

            # Soft noise reduction: remove small-magnitude noise
            threshold = 0.01 * np.max(np.abs(f_hat_filtered))
            f_hat_filtered[np.abs(f_hat_filtered) < threshold] = 0

            grad_A[..., i, j] = np.real(ifft2(f_hat_filtered))

    dot = np.sum(grad_u * grad_A)
    norm_grad_u = np.linalg.norm(grad_u)
    norm_A = np.linalg.norm(grad_A)
    Q = dot / (norm_grad_u * norm_A + 1e-10) if norm_grad_u > 1e-10 and norm_A > 1e-10 else 1.0
    KE = 0.5 * np.sum(u_grid**2 + v_grid**2)

    print(f"\n🧠 Q(t) Diagnostic (Cleaned) — {time_exact} — {level} hPa")
    print(f"📍 Region: {lat_min:.2f}–{lat_max:.2f}°N, {lon_min:.2f}–{lon_max:.2f}°W | kc = {kc}")
    print(f"Q(t) = {Q:.6f}  | KE = {KE:.2f}")
    print(f"‖∇u‖ = {norm_grad_u:.2e},  ‖A‖ = {norm_A:.2e},  ⟨∇u, A⟩ = {dot:.2e}")


