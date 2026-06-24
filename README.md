## AA-Stacked Silicene: Tight-Binding Band Structure Under Electric Fields

#### Abstract

This project presents a tight-binding computational study of AA-stacked silicene, focusing on its electronic band structure under both uniform and spatially modulated electric fields. Using a 4×4 Hamiltonian including intralayer and interlayer hopping terms, we compute energy bands along high-symmetry paths in the Brillouin zone and visualize band dispersion and k-space energy landscapes.

---
#### 1. Introduction

Silicene, a two-dimensional allotrope of silicon with a honeycomb lattice, exhibits graphene-like electronic properties with enhanced spin-orbit coupling and tunable band structure. In this work, we extend a tight-binding model to an AA-stacked bilayer configuration and investigate the effect of external electric fields on its electronic spectrum.

---
#### 2. Model Description

The system is modeled as a honeycomb lattice with lattice constant a = 2.28 Å. We consider AA stacking, where each atom in layer 1 is aligned vertically with its counterpart in layer 2. The Tight-Binding Hamiltonian includes:

- Intralayer hopping (t₀)
- Interlayer vertical hopping (t₁)
- Interlayer skew hopping (t₂)

The system is described by a 4×4 Hamiltonian in the basis: (A1, B1, A2, B2). The Bloch Hamiltonian depends on the structure factor: f(k) = Σ exp(i k · δᵢ) where δᵢ are nearest-neighbor vectors.

---
#### 3. Electric Field Effects

We consider two types of external perturbations: a uniform electric field (a constant on-site potential: V₀ applied equally to all sublattices), and a periodic electric field (a spatially modulated field: V(k) = V₀ cos(q · k), introducing k-dependent band modulation). 

---
#### 4. Computational Methods

Band structures are computed by:

1. Constructing H(k) on a 4×4 basis
2. Diagonalizing using numpy.linalg.eigvalsh
3. Sampling k-points along high-symmetry paths: Γ → M → K → Γ
4. Extending analysis to full k-space grids around K

---
#### 5. Results

The model produces four energy bands corresponding to bilayer sublattice degrees of freedom. We compute full 2D dispersion around the K point, revealing band splitting due to interlayer coupling and anisotropic curvature near Dirac points. The main effects due to the electric fields are:

- Uniform field shifts all bands rigidly
- Periodic field introduces k-dependent band modulation and symmetry breaking

---
#### 6. Key Findings

- Interlayer coupling (t₁, t₂) significantly modifies band splitting
- Electric fields introduce tunability in band structure
- AA stacking preserves high symmetry but lifts degeneracies

---
#### 7. Tools

- Python
- NumPy
- Matplotlib

---
#### 8. Conclusion
This work provides a computational framework for studying AA-stacked silicene under external fields using a tight-binding approach. The results highlight the tunability of electronic structure via interlayer coupling and external perturbations, relevant for 2D material engineering.

---
#### Acknowledgements

Undergraduate Computational Physics coursework in condensed matter theory.
