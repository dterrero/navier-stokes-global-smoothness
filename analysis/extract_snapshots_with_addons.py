import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
import imageio
import vtk
from vtk.util import numpy_support

# === SETTINGS ===
h5_path = "convection_simulation.h5"
output_dir = "extracted_snapshots"
snapshot_indices = [0, 499, 999]  # You can change
z_slice = 64  # For 2D slice plots
make_gif = True
make_vtk = True

os.makedirs(output_dir, exist_ok=True)
gif_frames = []

# === RAW VTK SAVE FUNCTION ===
def save_vti_numpy_fields(u, omega, u_mag, omega_mag, dx, filename):
    N = u.shape[0]
    image = vtk.vtkImageData()
    image.SetDimensions(N, N, N)
    image.SetSpacing(dx, dx, dx)
    image.SetOrigin(0, 0, 0)

    def add_field(name, array):
        flat_array = array.reshape(-1, array.shape[-1] if array.ndim == 4 else 1)
        vtk_array = numpy_support.numpy_to_vtk(flat_array, deep=True)
        vtk_array.SetName(name)
        image.GetPointData().AddArray(vtk_array)  # ✅ Use point data instead of cell data

    add_field("velocity", u)
    add_field("vorticity", omega)
    add_field("|u|", u_mag)
    add_field("|omega|", omega_mag)

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(image)
    writer.Write()


# === MAIN LOOP ===
with h5py.File(h5_path, "r") as f:
    velocity = f["velocity"]
    time_array = f["time"][:] if "time" in f else np.arange(len(velocity)) * 1e-4

    for i in snapshot_indices:
        print(f"\n📦 Processing snapshot {i}...")

        u = velocity[i]
        t = time_array[i]
        N = u.shape[0]
        dx = 2 * np.pi / N

        # === Compute vorticity ===
        ux, uy, uz = u[..., 0], u[..., 1], u[..., 2]
        curl_x = np.gradient(uz, dx, axis=1) - np.gradient(uy, dx, axis=2)
        curl_y = np.gradient(ux, dx, axis=2) - np.gradient(uz, dx, axis=0)
        curl_z = np.gradient(uy, dx, axis=0) - np.gradient(ux, dx, axis=1)
        omega = np.stack([curl_x, curl_y, curl_z], axis=-1)

        u_mag = np.linalg.norm(u, axis=-1)
        omega_mag = np.linalg.norm(omega, axis=-1)

        # === Save raw arrays ===
        np.save(f"{output_dir}/velocity_snapshot_{i}.npy", u)
        np.save(f"{output_dir}/vorticity_snapshot_{i}.npy", omega)

        # === Visualization (2D slice) ===
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        im1 = axes[0].imshow(u_mag[:, :, z_slice], cmap='viridis', origin='lower')
        axes[0].set_title(f"Velocity Magnitude\nz={z_slice}, t={t:.2f}")
        plt.colorbar(im1, ax=axes[0])

        im2 = axes[1].imshow(omega_mag[:, :, z_slice], cmap='inferno', origin='lower')
        axes[1].set_title(f"Vorticity Magnitude\nz={z_slice}, t={t:.2f}")
        plt.colorbar(im2, ax=axes[1])

        plot_path = f"{output_dir}/viz_snapshot_{i}.png"
        plt.tight_layout()
        plt.savefig(plot_path)
        plt.close()
        print(f"✅ Saved PNG: {plot_path}")

        if make_gif:
            gif_frames.append(imageio.v2.imread(plot_path))

        # === Save VTK .vti (ParaView) ===
        if make_vtk:
            vtk_path = f"{output_dir}/snapshot_{i}.vti"
            save_vti_numpy_fields(u, omega, u_mag, omega_mag, dx, vtk_path)
            print(f"✅ Saved VTK: {vtk_path}")

# === GIF Animation ===
if make_gif and len(gif_frames) > 1:
    gif_path = os.path.join(output_dir, "velocity_vorticity_evolution.gif")
    imageio.mimsave(gif_path, gif_frames, duration=1.5)
    print(f"\n🎞️ Saved GIF animation: {gif_path}")
