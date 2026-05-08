# DG Methods for Advection-Diffusion in DOLFINx

[![Made with MyST](https://img.shields.io/badge/made%20with-myst-orange)](https://curve.space/examples/pixels)
[![DOI](https://zenodo.org/badge/1211554310.svg)](https://doi.org/10.5281/zenodo.19615612)

Reference notes for implementing a discontinuous Galerkin (DG) advection-diffusion solver in [FESTIM](https://github.com/festim-dev/FESTIM) using [DOLFINx](https://github.com/FEniCS/dolfinx). The document covers the DG weak formulation (SIPG for diffusion, upwind flux for advection), boundary condition treatment, and verification via the Method of Manufactured Solutions (MMS).

## Contents

- [**paper.md**](paper.md) — derivation of the DG weak form, discussion of CG stabilisation methods, and the complete formulation with boundary conditions
- [**notebooks/example_1/mwe_1.ipynb**](notebooks/example_1/mwe_1.ipynb) — MMS verification on the unit square: step-by-step assembly of the weak form in DOLFINx, solution comparison against the exact solution, and mesh/degree convergence tests
- [**notebooks/example_2/mwe_2.ipynb**](notebooks/example_2/mwe_2.ipynb) — pipe flow with a half-Poiseuille velocity profile and a uniform volumetric source: boundary condition tagging, mass balance verification via consistent fluxes, and a Péclet number sweep showing the transition from diffusion-dominated to advection-dominated transport
- [**notebooks/example_3/mwe_3.ipynb**](notebooks/example_3/mwe_3.ipynb) — Stokes flow in a box geometry meshed in SALOME: Taylor-Hood Stokes solve to obtain a divergence-free velocity field, checkpoint-based coupling to the advection-diffusion solver, and flux balance verification on a non-trivial flow path
- [**notebooks/example_4/mwe_4.ipynb**](notebooks/example_4/mwe_4.ipynb) — OpenFOAM velocity field in a 3D box: reading a converged incompressible Navier-Stokes solution via foam2dolfinx, writing the mesh and velocity to an io4dolfinx checkpoint, and solving the advection-diffusion equation on the same 3D unstructured mesh with a do-nothing outlet condition and flux balance verification

## Getting Started

```bash
conda env create -f environment.yml
conda activate adv-diff
myst start
```

## Built With

[MyST Markdown](https://mystmd.org) — the paper and notebooks are rendered as an interactive article site.
