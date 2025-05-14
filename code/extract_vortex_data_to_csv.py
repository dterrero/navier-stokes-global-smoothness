import h5py
import csv

# Path to your .h5 file
h5_file = 'vortex_data_128.h5'
csv_file = 'vortex_summary.csv'

with h5py.File(h5_file, 'r') as hf, open(csv_file, 'w', newline='') as out_csv:
    writer = csv.writer(out_csv)
    writer.writerow(['step', 'Q', 'epsilon', 'energy'])

    time_data = hf['time_data']
    for step_key in sorted(time_data.keys()):
        step_grp = time_data[step_key]
        Q = step_grp.attrs['Q']
        epsilon = step_grp.attrs['epsilon']
        energy = step_grp.attrs['energy']
        step = int(step_key.split('_')[-1])
        writer.writerow([step, Q, epsilon, energy])

print(f"✅ CSV written to {csv_file}")
