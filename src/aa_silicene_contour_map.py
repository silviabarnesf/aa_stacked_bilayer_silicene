import numpy as np
import matplotlib.pyplot as plt

# Silicene parameters
a = 2.28  # Bond length (Angstrom)
t0 = -1.6  # Intralayer hopping (A-B in the same layer)
t1 = -0.1  # Interlayer vertical hopping (A1-A2, B1-B2) for silicene
t2 = 0.15  # Interlayer skew hopping (A1-B2, A2-B1) for silicene

def silicene_hamiltonian_AA(kx, ky):
    # Silicene tight-binding function f(k)
    f_k = (np.exp(1j * kx * a) +
           np.exp(1j * (-kx * a / 2 + np.sqrt(3) * ky * a / 2)) +
           np.exp(1j * (-kx * a / 2 - np.sqrt(3) * ky * a / 2)))

    # Build the 4x4 Hamiltonian
    H = np.zeros((4, 4), dtype=complex)

    # Layer 1
    H[0, 1] = t0 * f_k  # A1-B1
    H[1, 0] = np.conj(H[0, 1])

    # Layer 2
    H[2, 3] = t0 * f_k  # A2-B2
    H[3, 2] = np.conj(H[2, 3])

    # Interlayer AA vertical hopping (A1-A2 and B1-B2)
    H[0, 2] = t1
    H[2, 0] = t1
    H[1, 3] = t1
    H[3, 1] = t1

    # Interlayer AA skew hopping (A1-B2 and B1-A2)
    H[3, 0] = t2 * f_k  # A1-B2
    H[2, 1] = t2 * f_k  # B1-A2
    H[1, 2] = np.conj(H[2, 1])  # A2-B1
    H[0, 3] = np.conj(H[3, 0])  # B2-A1

    return H

# kx and ky around K point
Kx = 2 * np.pi / (3 * a)
Ky = (2 * np.pi / (3 * a)) * 1 / np.sqrt(3)

delta_k = 0.5  # Zoom range around K point
kx_values = np.linspace(Kx - delta_k, Kx + delta_k, 100)
ky_values = np.linspace(Ky - delta_k, Ky + delta_k, 100)
kx, ky = np.meshgrid(kx_values, ky_values)

# Calculate eigenvalues
eigenvalues = np.array([
    np.linalg.eigvalsh(silicene_hamiltonian_AA(kx[i, j], ky[i, j]))
    for i in range(len(kx_values))
    for j in range(len(ky_values))
])

# Reshape eigenvalues for plotting
eigenvalues = eigenvalues.reshape(len(kx_values), len(ky_values), 4)

# Plotting Contour Map (example: uppermost band)
plt.figure(figsize=(8, 6))
contour_filled = plt.contourf(kx, ky, eigenvalues[:, :, 3], levels=10, cmap='inferno')

# Add colorbar and labels
plt.colorbar(contour_filled, label='Energy (eV)')
plt.title('AA-Stacked Silicene Band Structure Contour Map (Top Band)', fontsize=16)
plt.xlabel('$k_x$ (1/\\AA)', fontsize=14)
plt.ylabel('$k_y$ (1/\\AA)', fontsize=14)
plt.grid(True)
plt.show()
