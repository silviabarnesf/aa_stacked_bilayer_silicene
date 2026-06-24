# PERIODIC ELECTRIC FIELD on Silicene AA stacked

import numpy as np
import matplotlib.pyplot as plt

# Parameters
a = 2.28  # bond length
t0 = -1.6
t1 = -0.1
t2 = 0.15
Vo = 3  # Field amplitude (eV)

# Periodicity of the field
q_vector = np.array([2*np.pi/20, 2*np.pi/20])  # modulation wavevector

G = np.array([0, 0])
M = (2 * np.pi / (3 * a)) * np.array([1, 0])
K = (2 * np.pi / (3 * a)) * np.array([1, 1/np.sqrt(3)])
k_path = [G, M, K, G]
k_labels = ['$\Gamma$', 'M', 'K', '$\Gamma$']
n_kpoints = 500

def interpolate_kpath(path, n):
    kpts = []
    for i in range(len(path) - 1):
        for j in range(n):
            k = path[i] + (path[i+1] - path[i]) * j / n
            kpts.append(k)
    return np.array(kpts)

def f_k(k):
    delta = [
        np.array([a, 0]),
        np.array([-a/2, np.sqrt(3)*a/2]),
        np.array([-a/2, -np.sqrt(3)*a/2])
    ]
    return sum(np.exp(1j * k.dot(d)) for d in delta)

def periodic_potential(k):
    """Periodic electric field modulation based on k (proxy for real space)"""
    return Vo * np.cos(q_vector.dot(k))

def H_AA_periodic(k):
    fk = f_k(k)
    Vp = periodic_potential(k)  # varies with k

    H = np.zeros((4, 4), dtype=complex)

    H[0, 1] = t0 * fk
    H[1, 0] = np.conj(H[0,1])
    H[2, 3] = t0 * fk
    H[3, 2] = np.conj(H[2,3])

    H[0,2] = t1
    H[2,0] = t1
    H[1,3] = t1
    H[3,1] = t1

    H[3,0] = t2 * fk
    H[2,1] = t2 * fk
    H[1,2] = np.conj(H[2,1])
    H[0,3] = np.conj(H[3,0])

    H[0, 0] = Vp
    H[1, 1] = Vp
    H[2, 2] = -Vp
    H[3, 3] = -Vp

    return H

kpoints = interpolate_kpath(k_path, n_kpoints)
energies = []
for k in kpoints:
    energies.append(np.sort(np.linalg.eigvalsh(H_AA_periodic(k))))
energies = np.array(energies)

fig, ax = plt.subplots(figsize=(8, 6))
for i in range(4):
    ax.plot(energies[:, i], color='darkgreen', lw=2)

tick_pos = [int(pos) for pos in np.linspace(0, len(kpoints), len(k_labels))]
ax.set_xticks(tick_pos)
ax.set_xticklabels(k_labels, fontsize=20)
ax.set_ylabel('Energy (eV)', fontsize=20)
ax.set_title('Silicene AA - Periodic Electric Field', fontsize=20)
ax.axhline(0, color='gray', linestyle='--', lw=2)
ax.set_ylim(-8, 8)
ax.set_xlim(0, 3*n_kpoints)
ax.tick_params(axis='y', labelsize=16)
plt.tight_layout()
plt.show()
