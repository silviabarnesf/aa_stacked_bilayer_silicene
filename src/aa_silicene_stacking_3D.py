import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Silicene Lattice Constant
a = 2.28  # Silicene bond length in Angstroms
t0 = -1.6   # Intralayer hopping (A-B in the same layer)
t1 = -0.1   # Interlayer vertical hopping (A1-A2, B1-B2)
t2 = 0.15   # Interlayer skew hopping (A1-B2, A2-B1)

def silicene_hamiltonian(kx, ky):
    # Silicene tight-binding function f(k)
    fk = (np.exp(1j * kx * a) +
          np.exp(1j * (-kx * a / 2 + np.sqrt(3) * ky * a / 2)) +
          np.exp(1j * (-kx * a / 2 - np.sqrt(3) * ky * a / 2)))

    # Hamiltonian
    H = np.zeros((4, 4), dtype=complex)

    # Layer 1
    H[0, 1] = t0 * fk  # A1-B1
    H[1, 0] = np.conj(H[0, 1])

    # Layer 2
    H[2, 3] = t0 * fk  # A2-B2
    H[3, 2] = np.conj(H[2, 3])

    # Interlayer AA coupling (A1-A2 and B1-B2)
    H[0, 2] = t1
    H[2, 0] = t1
    H[1, 3] = t1
    H[3, 1] = t1

    # Interlayer AB coupling (A1-B2 and B1-A2)
    H[3, 0] = t2 * fk  # A1-B2
    H[2, 1] = t2 * fk  # B1-A2
    H[1, 2] = np.conj(H[2, 1])  # A2-B1
    H[0, 3] = np.conj(H[3, 0])  # B2-A1

    return H

# kx and ky around K point
Kx = 2 * np.pi / (3 * a)
Ky = (2 * np.pi / (3 * a)) * 1 / np.sqrt(3)
delta_k = 0.5  # zoom range around K point

kx_values = np.linspace(Kx - delta_k, Kx + delta_k, 200)
ky_values = np.linspace(Ky - delta_k, Ky + delta_k, 200)
kx, ky = np.meshgrid(kx_values, ky_values)

# Calculate eigenvalues
eigenvalues = np.array([np.linalg.eigvalsh(silicene_hamiltonian(kx[i, j], ky[i, j]))
                        for i in range(len(kx_values))
                        for j in range(len(ky_values))])

# Reshape eigenvalues into 2D grids for each band
eigenvalues = eigenvalues.reshape(len(kx_values), len(ky_values), 4)
band1 = eigenvalues[:, :, 0]
band2 = eigenvalues[:, :, 1]
band3 = eigenvalues[:, :, 2]
band4 = eigenvalues[:, :, 3]

# 3D Plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot each band
ax.plot_surface(kx, ky, band1, cmap='viridis', alpha=0.8)
ax.plot_surface(kx, ky, band2, cmap='viridis', alpha=0.8)
ax.plot_surface(kx, ky, band3, cmap='viridis', alpha=0.8)
ax.plot_surface(kx, ky, band4, cmap='viridis', alpha=0.8)

ax.set_title('Silicene Tight-Binding Band Structure (3D)', fontsize=16)
ax.set_xlabel('$k_x$', fontsize=14)
ax.set_ylabel('$k_y$', fontsize=14)
ax.set_zlabel('Energy (eV)', fontsize=14)

plt.tight_layout()
plt.show()

# Print apparent bandgap
print("Max of lower band:", np.max(band2))
print("Min of upper band:", np.min(band3))
print("Apparent band gap:", np.min(band3) - np.max(band2))
