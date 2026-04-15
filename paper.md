---
title: DG Methods for Advection-Diffusion in DOLFINx
short_title: lorem
description: lorem
license:
  content: CC-BY-SA-3.0
keywords:
  - keyword 1
  - keyword 2
  - keyword 3
export:
  - format: pdf
    template: lapreprint
---

# DG Methods for Advection-Diffusion in DOLFINx

This article documents the development of Discontinuous Galerkin (DG) methods for the steady-state advection-diffusion equation, implemented using [DOLFINx](https://fenicsproject.org/). The work is motivated by efforts to couple [OpenFOAM](https://www.openfoam.com/) velocity fields to [FESTIM](https://festim.readthedocs.io/) for hydrogen transport simulations, where advection-diffusion of tritium in a flowing coolant is of interest.

## Local Setup

Clone the repository and create the conda environment:

```bash
git clone https://github.com/jhdark/advection_diffusion
cd advection_diffusion
conda env create -f environment.yml
conda activate adv-diff
```

Build and preview locally:

```bash
myst start
```

Or build static HTML:

```bash
myst build --html
# output in _build/html/
```
