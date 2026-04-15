---
title: DG Methods for Advection-Diffusion in DOLFINx
short_title: DG Advection-Diffusion
description: Reference notes for implementing a DG-based advection-diffusion class in FESTIM, building from MMS verification through to OpenFOAM-coupled simulations.
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

## The Advection-Diffusion Equation

We consider the steady-state advection-diffusion equation for a scalar quantity $u$ (e.g. concentration or temperature):

$$
\nabla \cdot (\mathbf{w} u) - \nabla \cdot (D \nabla u) = f \quad \text{in } \Omega
$$

where:
- $\mathbf{w}$ is the advective velocity field (prescribed, e.g. from OpenFOAM)
- $D > 0$ is the diffusion coefficient
- $f$ is a source term
- $\Omega \subset \mathbb{R}^d$ is the computational domain

Boundary conditions on $\partial\Omega$ split into inflow ($\mathbf{w} \cdot \mathbf{n} < 0$) and outflow ($\mathbf{w} \cdot \mathbf{n} \geq 0$) parts, with Dirichlet data prescribed on the inflow.

## The Péclet Number

The local cell Péclet number characterises the relative strength of advection to diffusion:

$$
\text{Pe}_h = \frac{|\mathbf{w}| h}{2D}
$$

where $h$ is a local mesh size. It governs which physical mechanism dominates:

| Regime | $\text{Pe}_h$ | Character |
|---|---|---|
| Diffusion-dominated | $\ll 1$ | Smooth, elliptic |
| Mixed | $\sim 1$ | Transition |
| Advection-dominated | $\gg 1$ | Transport-like, potential for spurious oscillations |

The appropriate numerical formulation depends strongly on $\text{Pe}_h$.

## Discontinuous Galerkin Methods

Standard continuous Galerkin (CG) methods can produce spurious oscillations when $\text{Pe}_h \gg 1$. Discontinuous Galerkin (DG) methods offer several advantages:

- **Upwinding** is natural for the advection term via numerical flux choices
- **Local conservation** is satisfied element-wise
- **High-order accuracy** on unstructured meshes
- **Flexibility** in handling complex boundary conditions weakly

The price is a larger global stiffness matrix (more degrees of freedom for the same mesh) and the need to choose penalty parameters carefully.

## Variational Formulation

Starting from the strong form, multiplying by a test function $v$ and integrating by parts over each element $K$, then summing over all elements, leads to the DG weak form. The key components are:

**Diffusion — Symmetric Interior Penalty (SIPG):**

$$
a_\text{diff}(u, v) = \sum_K \int_K D\, \nabla u \cdot \nabla v\, \mathrm{d}x
  - \sum_{F \in \mathcal{F}_I \cup \mathcal{F}_D} \int_F D\, \{\!\!\{ \nabla u \}\!\!\} \cdot \llbracket v \rrbracket\, \mathrm{d}s
  - \sum_{F \in \mathcal{F}_I \cup \mathcal{F}_D} \int_F D\, \llbracket u \rrbracket \cdot \{\!\!\{ \nabla v \}\!\!\}\, \mathrm{d}s
  + \sum_{F \in \mathcal{F}_I \cup \mathcal{F}_D} \int_F \frac{\alpha D}{h}\, \llbracket u \rrbracket \cdot \llbracket v \rrbracket\, \mathrm{d}s
$$

**Advection — upwind flux:**

$$
a_\text{adv}(u, v) = -\sum_K \int_K \mathbf{w} u \cdot \nabla v\, \mathrm{d}x
  + \sum_{F \in \mathcal{F}_I} \int_F \hat{f}_\text{up}(u)\, \llbracket v \rrbracket\, \mathrm{d}s
  + \text{boundary flux terms}
$$

where $\hat{f}_\text{up}$ is the upwind numerical flux: it selects $u$ from the upwind side of each interior face.

The penalty parameter $\alpha = C \cdot p^2$ (with $p$ the polynomial degree and $C \sim 10$) must be large enough to ensure coercivity of the diffusion bilinear form.

## Further Reading

- Arnold et al. (2002), *Unified analysis of discontinuous Galerkin methods for elliptic problems*, SIAM J. Numer. Anal.
- Di Pietro & Ern (2012), *Mathematical Aspects of Discontinuous Galerkin Methods*, Springer
- [DOLFINx documentation](https://docs.fenicsproject.org/dolfinx/main/python/)
