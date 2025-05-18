
from herbie import Herbie
import xarray as xr

times = ["2021-12-10 17:00", "2021-12-10 18:00", "2021-12-10 19:00", "2021-12-10 20:00"]
levels = [925, 850, 700]
datasets = []

for t in times:
    print(f"📥 Downloading: {t}...")
    H = Herbie(t, model="hrrr", product="prs", fxx=0, save_dir="hrrr_download", verbose=False)
    try:
        ds = H.xarray("UGRD:925 mb|VGRD:925 mb|UGRD:850 mb|VGRD:850 mb|UGRD:700 mb|VGRD:700 mb")
        datasets.append(ds)
    except Exception as e:
        print(f"⚠️ Failed to load {t}: {e}")

if datasets:
    print("🧪 Combining all datasets...")
    ds_merged = xr.concat(datasets, dim="time")
    ds_merged.to_netcdf("hrrr_mayfield_full.nc")
    print("✅ Saved as hrrr_mayfield_full.nc")
else:
    print("❌ No datasets to merge.")
