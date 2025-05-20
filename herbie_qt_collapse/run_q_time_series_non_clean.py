# ------------------------------------------------------------------------------
# run_q_time_series_non_clean.py
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
from scipy.spatial import cKDTree
from scipy.interpolate import RegularGridInterpolator
from datetime import datetime

def run_q_diagnostic(ds, time_str, level, lat_min, lat_max, lon_min, lon_max, kc):
    ds_slice = ds.sel(time=np.datetime64(time_str), isobaricInhPa=level, method="nearest")

    u = ds_slice['u'].values
    v = ds_slice['v'].values
    lat = ds_slice['latitude'].values
    lon = ds_slice['longitude'].values

    npts = 128
    lat_new = np.linspace(lat_min, lat_max, npts)
    lon_new = np.linspace(lon_min, lon_max, npts)
    lat_grid, lon_grid = np.meshgrid(lat_new, lon_new, indexing='ij')

    latlon_pairs = np.column_stack((lat.ravel(), lon.ravel()))
    tree = cKDTree(latlon_pairs)
    query_points = np.column_stack((lat_grid.ravel(), lon_grid.ravel()))
    _, idx = tree.query(query_points)
    y_idx, x_idx = np.unravel_index(idx, lat.shape)
    interp_points = np.column_stack((y_idx, x_idx))

    u_interp = RegularGridInterpolator((np.arange(u.shape[0]), np.arange(u.shape[1])), u, bounds_error=False, fill_value=np.nan)
    v_interp = RegularGridInterpolator((np.arange(v.shape[0]), np.arange(v.shape[1])), v, bounds_error=False, fill_value=np.nan)
    u_grid = u_interp(interp_points).reshape((npts, npts))
    v_grid = v_interp(interp_points).reshape((npts, npts))

    vel = np.stack([u_grid, v_grid], axis=-1)
    dy_vals = np.abs(np.diff(lat_grid, axis=0))
    dx_vals = np.abs(np.diff(lon_grid, axis=1))
    dy = np.nanmean(dy_vals)
    dx = np.nanmean(dx_vals)

    ny, nx, _ = vel.shape
    grad_u = np.zeros((ny, nx, 2, 2))
    for i in range(2):
        for j, spacing in enumerate([dy, dx]):
            grad_u[..., i, j] = np.gradient(vel[..., i], spacing, axis=j)

    KX = fftfreq(nx, d=dx) * 2 * np.pi
    KY = fftfreq(ny, d=dy) * 2 * np.pi
    KX, KY = np.meshgrid(KX, KY)
    Kmag = np.sqrt(KX**2 + KY**2)
    filter_mask = (Kmag <= kc)

    grad_A = np.zeros_like(grad_u)
    for i in range(2):
        for j in range(2):
            f_hat = fft2(grad_u[..., i, j])
            f_hat_filtered = f_hat * filter_mask
            grad_A[..., i, j] = np.real(ifft2(f_hat_filtered))

    dot = np.sum(grad_u * grad_A)
    norm_grad_u = np.linalg.norm(grad_u)
    norm_A = np.linalg.norm(grad_A)
    Q = dot / (norm_grad_u * norm_A + 1e-10) if norm_grad_u > 1e-10 and norm_A > 1e-10 else 1.0
    KE = 0.5 * np.sum(u_grid**2 + v_grid**2)

    print(f"\n🧠 Q(t) Diagnostic — {time_str} — {level} hPa")
    print(f"📍 Region: {lat_min:.2f}–{lat_max:.2f}°N, {lon_min:.2f}–{lon_max:.2f}°W | kc = {kc}")
    print(f"Q(t) = {Q:.6f}  | KE = {KE:.2f}")
    print(f"‖∇u‖ = {norm_grad_u:.2e},  ‖A‖ = {norm_A:.2e},  ⟨∇u, A⟩ = {dot:.2e}")

# === Load full dataset with time and multiple pressure levels
ds = xr.open_dataset("hrrr_mayfield_full.nc")

# Normalize longitudes if needed
if ds.longitude.max() > 180:
    ds = ds.assign_coords(longitude=ds.longitude.where(ds.longitude <= 180, ds.longitude - 360))

# Mayfield tornado outbreak location
lat_min, lat_max = 36.6, 37.0
lon_min, lon_max = -89.5, -89.0
kc = 10

# Times and levels to analyze
times = ["2021-12-10T17:00", "2021-12-10T18:00", "2021-12-10T19:00", "2021-12-10T20:00"]
levels = [925, 850, 700]  # hPa

for t in times:
    for lev in levels:
        run_q_diagnostic(ds, t, lev, lat_min, lat_max, lon_min, lon_max, kc)
