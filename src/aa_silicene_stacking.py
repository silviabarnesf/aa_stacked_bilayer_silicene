import numpy as np
import matplotlib.pyplot as plt

# Silicene Lattice Constant (bond length)
a = 2.28  # in Angstroms
t = -1.6  # Nearest-neighbor hopping energy in eV
t1 = -0.1  # Interlayer A-A, B-B hopping energy
t2 = 0.15  # Interlayer A-B, B-A hopping energy

# High-symmetry points in the Brillouin zone (Gamma, K, M, Gamma)
G = np.array([0, 0])
M = (2 * np.pi / (3 * a)) * np.array([1, 0])
K = (2 * np.pi / (3 * a)) * np.array([1, 1 / np.sqrt(3)])

# Path through BZ: G-M-K-G
k_path = [G, M, K, G]
k_labels = ['$\Gamma$', 'M', 'K', '$\Gamma$']

def interpolate_kpath(path, n):
    kpts = []
    for i in range(len(path) - 1):
        for j in range(n):
            k = path[i] + (path[i + 1] - path[i]) * j / n
            kpts.append(k)
    return np.array(kpts)

n_kpoints = 500

# f(k) = nearest-neighbor phase sum
def f_k(k):
    delta = [
        np.array([a, 0]),
        np.array([-a / 2, np.sqrt(3) * a / 2]),
        np.array([-a / 2, -np.sqrt(3) * a / 2])
    ]
    return sum(np.exp(1j * k.dot(d)) for d in delta)

# Build 4x4 tight-binding Hamiltonian
def H_AA(k):
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

    # Interlayer AB coupling (A1-B2 and B1-A2)
    H[3, 0] = t2 * fk
    H[2, 1] = t2 * fk
    H[1, 2] = np.conj(H[2, 1])
    H[0, 3] = np.conj(H[3, 0])

    return H

# Compute band structure
kpoints = interpolate_kpath(k_path, n_kpoints)
energies = []

for k in kpoints:
    H = H_AA(k)
    E = np.linalg.eigvalsh(H)
    energies.append(np.sort(E.real))
energies = np.array(energies)

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
for i in range(4):
    ax.plot(energies[:, i], label=f'Band {i+1}', color='blue', lw=2)

tick_pos = [0]
for i in range(1, len(k_path)):
    dk = np.linalg.norm(k_path[i] - k_path[i - 1])
    tick_pos.append(tick_pos[-1] + dk * n_kpoints)

tick_pos = [int(pos) for pos in np.linspace(0, len(kpoints), len(k_labels))]
ax.set_xticks(tick_pos)
ax.set_xticklabels(k_labels, fontsize=20)
ax.set_ylabel('Energy (eV)', fontsize=20, labelpad=10)
ax.set_title('Band Structure of AA-stacked Silicene', fontsize=20)
ax.axhline(0, color='gray', linestyle='--', lw=2)
ax.set_ylim(-8, 8)
ax.set_xlim(0, 3 * n_kpoints)
ax.tick_params(axis='y', labelsize=16)

for spine in ax.spines.values():
    spine.set_linewidth(2)

plt.tight_layout()
plt.show()
