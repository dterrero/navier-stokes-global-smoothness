import h5py
import numpy as np
import pandas as pd

h5_path = "convection_simulation.h5"
csv_path = "convection_diagnostics_128.csv"

with h5py.File(h5_path, "r") as f:
    data = {}
    keys = ["time", "Nu", "heat_flux", "Q", "Q_num", "norm_grad_u", "norm_A", "KE"]
    for key in keys:
        if key in f:
            data[key] = f[key][:]
        else:
            print(f"Warning: '{key}' not found. Filling with NaNs.")
            data[key] = None

# Determine minimum length of all present arrays
lengths = [len(v) for v in data.values() if v is not None]
min_len = min(lengths)

# Truncate all arrays to min_len and fill missing with NaNs
for key in data:
    if data[key] is None:
        data[key] = np.full(min_len, np.nan)
    else:
        data[key] = data[key][:min_len]

# Create and save the DataFrame
df = pd.DataFrame(data)
df.to_csv(csv_path, index=False)
print(f"Diagnostics saved to {csv_path} (length = {min_len})")
