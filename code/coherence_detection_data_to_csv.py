
import h5py
import numpy as np
import pandas as pd

# Load HDF5 file
with h5py.File("coherence_detection.h5", "r") as f:
    time = f["time_full"][:]
    Q = f["Q"][:]
    Nu = f["Nu"][:]
    heat_flux = f["heat_flux"][:]

# Construct DataFrame
df = pd.DataFrame({
    "time": time,
    "Q": Q,
    "Nu": Nu,
    "heat_flux": heat_flux
})

# Save to CSV
df.to_csv("coherence_detection.csv", index=False)
print("CSV exported to 'coherence_detection.csv'")
