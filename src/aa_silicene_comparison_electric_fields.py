import numpy as np
import matplotlib.pyplot as plt

# AA Silicene Lattice Constants
a = 1.42  # Lattice constant (Å)
t = -2.7  # Nearest-neighbor hopping energy in eV
t1 = -0.263
t2 = 0.32
Vo = 5  # Stronger electric field (eV)

# High-symmetry points in the Brillouin zone (Gamma, K, M, Gamma)
G = np.array([0, 0])
M = (2 * np.pi / (3 * a)) * np.array([1, 0])  # middle point
K = (2 * np.pi / (3 * a)) * np.array([1, 1 / np.sqrt(3)])  # corner point

# Path through BZ: G-K-M-G
k_path = [G, M, K, G]
k_labels = ['$\Gamma$', 'M', 'K', '$\Gamma$']

# k-point interpolation
def interpolate_kpath(path, n):
    kpts = []
    for i in range(len(path) - 1):
        for j in range(n):
            k = path[i] + (path[i + 1] - path[i]) * j / n
            kpts.append(k)
    return np.array(kpts)

n_kpoints = 200

# f(k) = nearest-neighbor phase sum
def f_k(k):
    delta = [
        np.array([a, 0]),  # sigma_1
        np.array([-a / 2, np.sqrt(3) * a / 2]),  # sigma_2
        np.array([-a / 2, -np.sqrt(3) * a / 2])  # sigma_3
    ]
    return sum(np.exp(1j * k.dot(d)) for d in delta)

# Build 4x4 tight-binding Hamiltonian
def H_AA(k, V0=0, q_vector=None):
    fk = f_k(k)
    H = np.zeros((4, 4), dtype=complex)

    # Layer 1
    H[0, 1] = t * fk  # A1-B1
    H[1, 0] = np.conj(H[0, 1])

    # Layer 2
    H[2, 3] = t * fk  # A2-B2
    H[3, 2] = np.conj(H[2, 3])

    # Interlayer AA coupling (A1-A2 and B1-B2)
    H[0, 2] = t1
    H[2, 0] = t1
    H[1, 3] = t1
    H[3, 1] = t1

    # Interlayer AA coupling (A1-B2 and B1-A2)
    H[3, 0] = t2 * fk  # A1-B2
    H[2, 1] = t2 * fk  # B1-A2
    H[1, 2] = np.conj(H[2, 1])  # A2-B1
    H[0, 3] = np.conj(H[3, 0])  # B2-A1

    # Add electric field (uniform + periodic)
    if q_vector is not None:
        field = V0 * np.cos(k.dot(q_vector))
        H[0, 0] += field
        H[1, 1] += field
        H[2, 2] += field
        H[3, 3] += field
    else:
        # Uniform field
        H[0, 0] += V0
        H[1, 1] += V0
        H[2, 2] += V0
        H[3, 3] += V0

    return H

# Compute band structure for periodic and uniform electric field
kpoints = interpolate_kpath(k_path, n_kpoints)
energies_uniform = []
energies_periodic = []

# Periodic electric field (shorter q_vector for a stronger modulation)
q_vector = np.array([2*np.pi/5, 2*np.pi/5])  # Shorter periodicity (stronger oscillation)

for k in kpoints:
    H_uniform = H_AA(k, V0=5)  # Uniform electric field
    E_uniform = np.linalg.eigvalsh(H_uniform)
    energies_uniform.append(np.sort(E_uniform.real))

    H_periodic = H_AA(k, V0=5, q_vector=q_vector)  # Periodic electric field
    E_periodic = np.linalg.eigvalsh(H_periodic)
    energies_periodic.append(np.sort(E_periodic.real))

energies_uniform = np.array(energies_uniform)
energies_periodic = np.array(energies_periodic)

# Plotting
fig, ax = plt.subplots(figsize=(8, 6))
for i in range(4):
    ax.plot(energies_uniform[:, i], label=f'Band {i+1} (Uniform)', color='blue', lw=2)
for i in range(4):
    ax.plot(energies_periodic[:, i], label=f'Band {i+1} (Periodic)', color='red', linestyle='--', lw=2)

# Add symmetry point markers
tick_pos = [0]
for i in range(1, len(k_path)):
    dk = np.linalg.norm(k_path[i] - k_path[i - 1])
    tick_pos.append(tick_pos[-1] + dk * n_kpoints)

tick_pos = [int(pos) for pos in np.linspace(0, len(kpoints), len(k_labels))]
ax.set_xticks(tick_pos)
ax.set_xticklabels(k_labels, fontsize=16)
ax.set_ylabel('Energy (eV)', fontsize=16)
ax.set_title('Band Structure of Silicene with Electric Fields', fontsize=18)
ax.axhline(0, color='gray', linestyle='--', lw=2)  # Fermi level
ax.set_ylim(-8, 8)
ax.set_xlim(0, 3 * n_kpoints)
ax.tick_params(axis='y', labelsize=14)

# Clean up plot aesthetics
for spine in ax.spines.values():
    spine.set_linewidth(2)

plt.tight_layout()
plt.legend()
plt.show()
