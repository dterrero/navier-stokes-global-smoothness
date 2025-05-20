# ------------------------------------------------------------------------------
# run_q_time_series_original_vs_artifacts_removal.py
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
import xarray as xr
from numpy.fft import fft2, ifft2, fftfreq
from scipy.interpolate import RegularGridInterpolator
import pandas as pd

results = []

def run_q_diagnostic(ds, time_exact, level, lat_min, lat_max, lon_min, lon_max, kc):
    try:
        # Select exact time match only — avoid duplication
        ds_slice = ds.sel(time=np.datetime64(time_exact))
    except KeyError:
        print(f"⚠️ Time {time_exact} not found in dataset.")
        return

    # Filter by level
    ds_slice = ds_slice.sel(isobaricInhPa=level, method="nearest")

    u = ds_slice['u'].values
    v = ds_slice['v'].values
    lat = ds_slice['latitude'].values
    lon = ds_slice['longitude'].values

    npts = 64
    lat_new = np.linspace(lat_min, lat_max, npts)
    lon_new = np.linspace(lon_min, lon_max, npts)
    lat_grid, lon_grid = np.meshgrid(lat_new, lon_new, indexing='ij')
    points = np.column_stack((lat_grid.ravel(), lon_grid.ravel()))

    # Interpolation
    u_interp = RegularGridInterpolator((lat[:, 0], lon[0, :]), u, bounds_error=False, fill_value=np.nan)
    v_interp = RegularGridInterpolator((lat[:, 0], lon[0, :]), v, bounds_error=False, fill_value=np.nan)

    u_grid = u_interp(points).reshape((npts, npts))
    v_grid = v_interp(points).reshape((npts, npts))

    vel = np.stack([u_grid, v_grid], axis=-1)
    dy = np.nanmean(np.abs(np.diff(lat_grid, axis=0)))
    dx = np.nanmean(np.abs(np.diff(lon_grid, axis=1)))

    grad_u = np.zeros((npts, npts, 2, 2))
    for i in range(2):
        for j, spacing in enumerate([dy, dx]):
            grad_u[..., i, j] = np.gradient(vel[..., i], spacing, axis=j)

    KX = fftfreq(npts, d=dx) * 2 * np.pi
    KY = fftfreq(npts, d=dy) * 2 * np.pi
    KX, KY = np.meshgrid(KX, KY)
    Kmag = np.sqrt(KX**2 + KY**2)
    filter_mask = (Kmag <= kc)

    grad_A = np.zeros_like(grad_u)
    for i in range(2):
        for j in range(2):
            f_hat = fft2(grad_u[..., i, j])
            grad_A[..., i, j] = np.real(ifft2(f_hat * filter_mask))

    dot = np.sum(grad_u * grad_A)
    norm_grad_u = np.linalg.norm(grad_u)
    norm_A = np.linalg.norm(grad_A)
    Q = dot / (norm_grad_u * norm_A + 1e-10) if norm_grad_u > 1e-10 and norm_A > 1e-10 else 1.0
    KE = 0.5 * np.sum(u_grid**2 + v_grid**2)

    print(f"\n🧠 Q(t) Diagnostic — {time_exact} — {level} hPa")
    print(f"📍 Region: {lat_min:.2f}–{lat_max:.2f}°N, {lon_min:.2f}–{lon_max:.2f}°W | kc = {kc}")
    print(f"Q(t) = {Q:.6f}  | KE = {KE:.2f}")
    print(f"‖∇u‖ = {norm_grad_u:.2e},  ‖A‖ = {norm_A:.2e},  ⟨∇u, A⟩ = {dot:.2e}")
    return Q  # ✅ This is the missing piece

def run_q_diagnostic_cleaned(ds, time_exact, level, lat_min, lat_max, lon_min, lon_max, kc):
    try:
        ds_slice = ds.sel(time=np.datetime64(time_exact), isobaricInhPa=level, method="nearest")
    except KeyError:
        print(f"⚠️ Time {time_exact} not found in dataset.")
        return

    u = ds_slice['u'].values
    v = ds_slice['v'].values
    lat = ds_slice['latitude'].values
    lon = ds_slice['longitude'].values

    npts = 64
    pad = 4
    npts_total = npts + 2 * pad
    lat_full = np.linspace(lat_min, lat_max, npts_total)
    lon_full = np.linspace(lon_min, lon_max, npts_total)
    lat_grid, lon_grid = np.meshgrid(lat_full, lon_full, indexing='ij')
    points = np.column_stack((lat_grid.ravel(), lon_grid.ravel()))

    u_interp = RegularGridInterpolator((lat[:, 0], lon[0, :]), u, bounds_error=False, fill_value=np.nan)
    v_interp = RegularGridInterpolator((lat[:, 0], lon[0, :]), v, bounds_error=False, fill_value=np.nan)
    u_grid_full = u_interp(points).reshape((npts_total, npts_total))
    v_grid_full = v_interp(points).reshape((npts_total, npts_total))

    u_grid = u_grid_full[pad:-pad, pad:-pad]
    v_grid = v_grid_full[pad:-pad, pad:-pad]
    vel = np.stack([u_grid, v_grid], axis=-1)

    dy = np.nanmean(np.abs(np.diff(lat_full[pad:-pad])))
    dx = np.nanmean(np.abs(np.diff(lon_full[pad:-pad])))

    grad_u = np.zeros((npts, npts, 2, 2))
    for i in range(2):
        for j, spacing in enumerate([dy, dx]):
            grad_u[..., i, j] = np.gradient(vel[..., i], spacing, axis=j)

    KX = fftfreq(npts, d=dx) * 2 * np.pi
    KY = fftfreq(npts, d=dy) * 2 * np.pi
    KX, KY = np.meshgrid(KX, KY)
    Kmag = np.sqrt(KX**2 + KY**2)
    mask = Kmag <= kc

    grad_A = np.zeros_like(grad_u)
    for i in range(2):
        for j in range(2):
            F = fft2(grad_u[..., i, j])
            threshold = 0.01 * np.max(np.abs(F))
            F[np.abs(F) < threshold] = 0  # Remove filament noise
            grad_A[..., i, j] = np.real(ifft2(F * mask))

    dot = np.sum(grad_u * grad_A)
    norm_grad_u = np.linalg.norm(grad_u)
    norm_A = np.linalg.norm(grad_A)
    Q = dot / (norm_grad_u * norm_A + 1e-10)
    print(f"🧹 Cleaned Q(t) = {Q:.6f}")
    return Q

# === Load dataset
ds = xr.open_dataset("hrrr_mayfield_full.nc")

# Normalize longitudes
if ds.longitude.max() > 180:
    ds = ds.assign_coords(longitude=ds.longitude.where(ds.longitude <= 180, ds.longitude - 360))

print("\n🕒 Available time values in dataset:")
print(ds.time.values)

# Parameters
kc = 10
lat_min, lat_max = 36.6, 37.0
lon_min, lon_max = -89.5, -89.0

times = ["2021-12-10T17:00", "2021-12-10T18:00", "2021-12-10T19:00", "2021-12-10T20:00"]
levels = [925, 850, 700]
kc = 10
lat_min, lat_max = 36.6, 37.0
lon_min, lon_max = -89.5, -89.0

for time_str in times:
    for level in levels:
        print(f"\n=== Processing {time_str} | {level} hPa ===")

        try:
            q_orig = run_q_diagnostic(ds, time_str, level, lat_min, lat_max, lon_min, lon_max, kc)
            q_clean = run_q_diagnostic_cleaned(ds, time_str, level, lat_min, lat_max, lon_min, lon_max, kc)
        except Exception as e:
            print(f"⚠️ Failed for {time_str}, {level}: {e}")
            continue

        if q_orig is not None and q_clean is not None and not np.isnan(q_orig) and not np.isnan(q_clean):
            results.append({
                "Time": time_str,
                "Level_hPa": level,
                "Q_Original": round(q_orig, 6),
                "Q_Cleaned": round(q_clean, 6),
                "Delta_Q": round(q_clean - q_orig, 6)
            })
        else:
            print(f"⚠️ Skipping entry for {time_str}, {level} hPa due to missing Q value.")


df = pd.DataFrame(results)
csv_path = "q_diagnostics_comparison_mayfield.csv"
df.to_csv(csv_path, index=False)
print(f"\n✅ Results saved to: {csv_path}")

