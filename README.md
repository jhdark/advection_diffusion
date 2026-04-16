# DG Methods for Advection-Diffusion in DOLFINx

[![Made with MyST](https://img.shields.io/badge/made%20with-myst-orange)](https://curve.space/examples/pixels)
[![DOI](https://zenodo.org/badge/1211554310.svg)](https://doi.org/10.5281/zenodo.19615612)

Reference notes for implementing a discontinuous Galerkin (DG) advection-diffusion solver in [FESTIM](https://github.com/festim-dev/FESTIM) using [DOLFINx](https://github.com/FEniCS/dolfinx). The document covers the DG weak formulation (SIPG for diffusion, upwind flux for advection), boundary condition treatment, and verification via the Method of Manufactured Solutions (MMS).

## Contents

- [**paper.md**](paper.md) — derivation of the DG weak form, discussion of CG stabilisation methods, and the complete formulation with boundary conditions
- [**notebooks/mwe_1.ipynb**](notebooks/mwe_1.ipynb) — MMS verification on the unit square: step-by-step assembly of the weak form in DOLFINx, solution comparison against the exact solution, and mesh/degree convergence tests
- [**notebooks/mwe_2.ipynb**](notebooks/mwe_2.ipynb) — further examples

## Getting Started

```bash
conda env create -f environment.yml
conda activate adv-diff
myst start
```

## Built With

[MyST Markdown](https://mystmd.org) — the paper and notebooks are rendered as an interactive article site.
