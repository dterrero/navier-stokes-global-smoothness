# ------------------------------------------------------------------------------
# forced_shear_layer_with_coherence_disruption.py
#
# Computes the Coherence Quotient Q(t), Kinetic Energy (KE), and total
# gradient energy (proxy for dissipation) over time in a 3D shear layer.
# A coherent disruption (burst) is injected at t = 400 to analyze impact.
#
# Author: Dickson A. Terrero
# License: CC BY-NC 4.0
# https://creativecommons.org/licenses/by-nc/4.0/
# ------------------------------------------------------------------------------
# ⚠️ Usage Notice:
# This script and associated formulas are shared for **educational and research purposes only**.
# Commercial use is **not permitted** under the terms of the license.
#
# If you wish to use this method in a commercial product or service,
# please contact the author to discuss licensing terms.
# ------------------------------------------------------------------------------

# Parameters
N = 64
L = 2 * np.pi
dx = L / N
dt = 0.01
steps = 1000
kc = 15
viscosity = 1e-3

# Grids
x = np.linspace(0, L, N, endpoint=False)
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')

# FFT wave numbers
K = fftfreq(N, d=dx) * 2 * np.pi
KX, KY, KZ = np.meshgrid(K, K, K, indexing='ij')
K2 = KX**2 + KY**2 + KZ**2
K2[0, 0, 0] = 1e-10

# Initial velocity field (shear layer in x-direction)
u = np.zeros((3, N, N, N))
u[0] = np.tanh((Y - L / 2) * 10)

# Coherence Quotient Function
def compute_Q(grad_u, kc, return_diagnostics=False):
    grad_A = np.zeros_like(grad_u)
    K_mag = np.sqrt(KX**2 + KY**2 + KZ**2)
    filter_mask = (K_mag <= kc)

    for i in range(3):
        for j in range(3):
            grad_ij = grad_u[..., i, j]
            grad_hat = fftn(grad_ij)
            grad_hat_filtered = grad_hat * filter_mask
            grad_A[..., i, j] = np.real(ifftn(grad_hat_filtered))

    dot = np.sum(grad_u * grad_A)
    norm_u = np.linalg.norm(grad_u)
    norm_A = np.linalg.norm(grad_A)

    if norm_u < 1e-10 or norm_A < 1e-10:
        Q_val = 1.0
    else:
        Q_val = dot / (norm_u * norm_A + 1e-10)

    if return_diagnostics:
        return Q_val, dot, norm_u, norm_A
    else:
        return Q_val

# HDF5 output
with h5py.File("shear_layer_q_results.h5", "w") as h5f:
    dset_q = h5f.create_dataset("Q", (steps,), dtype='f')
    dset_ke = h5f.create_dataset("KE", (steps,), dtype='f')
    dset_grad_energy = h5f.create_dataset("GradientEnergy", (steps,), dtype='f')

    # Main simulation loop
    for step in tqdm(range(steps), desc="Simulating Coherence Dynamics"):
        # Inject energy burst at step 400
        if step == 400:
            burst = np.exp(-((X - L/2)**2 + (Y - L/2)**2 + (Z - L/2)**2) / 0.01)
            u[1] += burst * 5.0

        # Compute velocity gradient
        grad_u = np.zeros((N, N, N, 3, 3))
        for i in range(3):
            for j in range(3):
                grad_u[..., i, j] = np.gradient(u[i], dx, axis=j)

        # Compute diagnostics
        Q = compute_Q(grad_u, kc)
        KE = 0.5 * np.sum(u**2)
        grad_energy = np.sum(grad_u**2)

        # Store
        dset_q[step] = Q
        dset_ke[step] = KE
        dset_grad_energy[step] = grad_energy

        # Print every 50 steps
        if step % 50 == 0:
            print(f"Step {step:4d} | Q(t) = {Q:.5f} | KE = {KE:.2f} | ∇u² = {grad_energy:.2e}")
