---
title: DG Methods for Advection-Diffusion in DOLFINx
short_title: DG Advection-Diffusion
description: Reference notes for implementing a DG-based advection-diffusion class in FESTIM, building from MMS verification through to OpenFOAM-coupled simulations.
license:
  content: CC-BY-SA-3.0
keywords:
  - discontinuous Galerkin
  - advection-diffusion
  - FESTIM
  - DOLFINx
  - finite element methods
export:
  - format: pdf
    template: lapreprint
---

## The Advection-Diffusion Equation

We consider the steady-state advection-diffusion equation for a scalar quantity $u$, in this case in m$^{-3}$:

$$
\nabla \cdot (\mathbf{w} u) - \nabla \cdot (D \nabla u) = f \quad \text{in } \Omega
$$

where $\mathbf{w}$ is a prescribed advective velocity field in m$\,$s$^{-1}$, $D$ is the diffusion coefficient in m$^{2}\,$s$^{-1}$, $f$ is a source term in m$^{2}\,$s$^{-1}$, and $\Omega \subset \mathbb{R}^d$ is the computational domain.

The relative importance of the two transport mechanisms is captured by the cell Péclet number:

$$
\text{Pe} = \frac{|\mathbf{w}| h}{2D}
$$

where $h$ is a local mesh size. When $\text{Pe} \ll 1$ diffusion dominates and the problem is essentially elliptic; when $\text{Pe} \gg 1$ advection dominates and the solution can develop sharp internal layers and boundary layers.

## Numerical Challenges and the Case for DG

Standard continuous Galerkin (CG) finite element methods work well for the pure diffusion limit, but suffer from spurious non-physical oscillations as $\text{Pe}$ grows. These oscillations arise because CG lacks any intrinsic upwinding mechanism, as the scheme treats upwind and downwind information symmetrically, even though advection is inherently directional.

Several stabilisation strategies have been developed to address this within the CG framework:

- **Isotropic artificial diffusion** adds a scalar diffusion term $D_\text{art} = \delta h |\mathbf{w}|$ uniformly in all directions, where $\delta$ is a dimensionless tuning parameter and $h$ is the local mesh size. It is simple to implement but *inconsistent*, as the added diffusion modifies the original problem, introducing cross-stream smearing and reducing accuracy, particularly near sharp layers [(COMSOL, 2020)](https://www.comsol.com/blogs/understanding-stabilization-methods).
- **Streamline-upwind Petrov–Galerkin (SUPG)** and **Galerkin least-squares (GLS)** are *consistent* stabilisations: they add residual-weighted terms that introduce numerical diffusion strictly along the streamline direction. Because the added terms vanish when the exact solution is substituted, the convergence order is preserved. These are generally preferred over isotropic diffusion.

These methods can be effective, but all require tuning parameters, add complexity to the formulation, and are fundamentally workarounds for a framework not designed with advection-dominated transport in mind.

Discontinuous Galerkin (DG) methods address the root cause rather than patching the symptom. Because the approximation space allows inter-element discontinuities, information can be passed between elements through numerical fluxes, and the choice of flux naturally encodes upwinding for the advection term. This gives DG several structural advantages:

- **Natural upwinding**: the upwind flux for advection is a direct consequence of the DG framework, not an add-on stabilisation.
- **Local conservation**: the flux balance is satisfied element-wise, which matters for transport quantities.
- **Flexibility in boundary conditions**: all conditions are imposed weakly through face integrals, giving a uniform treatment.
- **High-order accuracy**: high-degree polynomial spaces can be used on unstructured meshes without additional complications.

The trade-off is a larger global system (DG has more degrees of freedom than CG for the same mesh) and the need to select a penalty parameter for the diffusion term carefully. Both are manageable in practice.

## Derivation of the DG Weak Form

### Notation

Let $\mathcal{T}_h = \{K\}$ be a triangulation of $\Omega$, and let $\mathcal{F}_I$ and $\mathcal{F}_B$ denote the sets of interior and boundary faces. The DG trial and test space is:

$$
V_h = \{v \in L^2(\Omega) : v|_K \in \mathbb{P}_p(K),\ \forall K \in \mathcal{T}_h\}
$$

For an interior face $F \in \mathcal{F}_I$ shared by elements $K_0$ and $K_1$, let $\mathbf{n} = \mathbf{n}_0$ be the outward unit normal from $K_0$. We define the average and jump operators for a scalar $q$ and vector $\boldsymbol{\tau}$ as:

$$
\langle q \rangle := \frac{q_0 + q_1}{2}, \qquad \llbracket q \rrbracket := q_0 - q_1
$$

$$
\langle \boldsymbol{\tau} \rangle := \frac{\boldsymbol{\tau}_0 + \boldsymbol{\tau}_1}{2}, \qquad \llbracket \boldsymbol{\tau} \rrbracket := \boldsymbol{\tau}_0 \cdot \mathbf{n}_0 + \boldsymbol{\tau}_1 \cdot \mathbf{n}_1
$$

Note that $\llbracket \boldsymbol{\tau} \rrbracket = (\boldsymbol{\tau}_0 - \boldsymbol{\tau}_1) \cdot \mathbf{n}$ is the scalar normal jump of a vector, while $\llbracket q \rrbracket = q_0 - q_1$ is the scalar jump of $q$ oriented with $\mathbf{n}$. On a boundary face there is a single element, so $\langle q \rangle = q$ and $\llbracket q \rrbracket = q$.

### Element-wise Integration by Parts

Multiply the strong form by $v \in V_h$ and integrate over a single element $K$:

$$
\int_K \nabla \cdot (\mathbf{w} u)\, v\, \mathrm{d}x - \int_K \nabla \cdot (D\nabla u)\, v\, \mathrm{d}x = \int_K f\, v\, \mathrm{d}x
$$

Applying the divergence theorem to each term separately:

$$
\begin{split}
&-\int_K (\mathbf{w} u) \cdot \nabla v\, \mathrm{d}x + \int_{\partial K} (\mathbf{w} \cdot \mathbf{n}_K)\, u\, v\, \mathrm{d}s \\
&+ \int_K D\nabla u \cdot \nabla v\, \mathrm{d}x - \int_{\partial K} D(\nabla u \cdot \mathbf{n}_K)\, v\, \mathrm{d}s = \int_K f\, v\, \mathrm{d}x
\end{split}
$$

Summing over all elements $K \in \mathcal{T}_h$ and assembling the $\partial K$ integrals across shared faces using the average and jump operators yields the global weak form.

### Diffusion Term: SIPG

The assembly identity for the diffusion boundary integrals gives:

$$
\begin{align}
    \sum_K \int_{\partial K} D(\nabla u \cdot \mathbf{n}_K)\, v\, \mathrm{d}s
    = &\sum_{F \in \mathcal{F}_I} \int_F \left( \langle D\nabla u \rangle \cdot \mathbf{n}\, \llbracket v \rrbracket
    + \llbracket D\nabla u \rrbracket\, \langle v \rangle \right) \mathrm{d}s \\
    &+ \sum_{F \in \mathcal{F}_B} \int_F D(\nabla u \cdot \mathbf{n})\, v\, \mathrm{d}s
\end{align}
$$

In the continuous setting the normal flux of a smooth solution is single-valued, so $\llbracket D\nabla u \rrbracket = 0$ and the second interior face term can be dropped from the bilinear form without compromising consistency, since it vanishes exactly when the true solution is substituted. In the discrete setting $u_h$ is genuinely discontinuous and $\llbracket D\nabla u_h \rrbracket$ is not zero in general, so omitting this term does modify the discrete system; stability is recovered through the penalty contribution introduced below, rather than by any smoothness of $u_h$. Subtracting the remaining single-valued flux term from the element-wise bulk integrals gives a consistent but non-symmetric form. SIPG adds two further face contributions to restore symmetry and coercivity:

1. **Symmetry**: $-\displaystyle\sum_{F \in \mathcal{F}_I} \int_F \langle D\nabla v \rangle \cdot \mathbf{n}\, \llbracket u \rrbracket\, \mathrm{d}s$, which is zero for the exact solution (consistency preserved) and symmetrises the bilinear form.
2. **Penalty**: $+\displaystyle\sum_{F \in \mathcal{F}_I} \int_F \dfrac{\alpha D}{h}\, \llbracket u \rrbracket \llbracket v \rrbracket\, \mathrm{d}s$, which penalises inter-element jumps and restores coercivity.

The interior SIPG bilinear form is therefore:

$$
\begin{align}
a_\text{diff}^\text{int}(u, v)
&= \sum_K \int_K D\, \nabla u \cdot \nabla v\, \mathrm{d}x \\
&\quad - \sum_{F \in \mathcal{F}_I} \int_F \langle D\nabla u \rangle \cdot \mathbf{n}\, \llbracket v \rrbracket\, \mathrm{d}s
   - \sum_{F \in \mathcal{F}_I} \int_F \langle D\nabla v \rangle \cdot \mathbf{n}\, \llbracket u \rrbracket\, \mathrm{d}s \\
&\quad + \sum_{F \in \mathcal{F}_I} \int_F \frac{\alpha D}{h}\, \llbracket u \rrbracket \llbracket v \rrbracket\, \mathrm{d}s
\end{align}
$$

The penalty parameter $\alpha = Cp^2$ (with $C \sim 10$ and $p$ the polynomial degree) must exceed a threshold that depends on the element geometry.

### Advection Term: Upwind Flux

After summing over elements, the advection face integrals over an interior face $F$ (with normal $\mathbf{n} = \mathbf{n}_0$) yield:

$$
\int_F (\mathbf{w} \cdot \mathbf{n})\bigl(u_0\, v_0 - u_1\, v_1\bigr)\, \mathrm{d}s
$$

Since $u$ is multi-valued on $F$, we replace it with the upwind value:

$$
u_\text{up} = \begin{cases} u_0 & \text{if } \mathbf{w} \cdot \mathbf{n} > 0 \\ u_1 & \text{if } \mathbf{w} \cdot \mathbf{n} \leq 0 \end{cases}
$$

which can be written compactly as:

$$
(\mathbf{w} \cdot \mathbf{n})\, u_\text{up} = \langle \mathbf{w} u \rangle \cdot \mathbf{n} + \tfrac{1}{2}|\mathbf{w} \cdot \mathbf{n}|\, \llbracket u \rrbracket
$$

The face integral then becomes $(\mathbf{w}\cdot\mathbf{n})\, u_\text{up}\, \llbracket v \rrbracket$, where $\llbracket v \rrbracket = v_0 - v_1$ is the oriented scalar jump of the test function.

The compact identity above assumes that $\mathbf{w}\cdot\mathbf{n}$ is single-valued on $F$, i.e. $\llbracket \mathbf{w}\cdot\mathbf{n} \rrbracket = 0$. This holds whenever $\mathbf{w}$ is prescribed analytically or is represented in a continuous finite element space, but it fails in general for a velocity field interpolated from an external solver into a DG space. The practical consequences, and how to handle them, are discussed in the section on coupling with an external velocity field.

The interior advection contribution is:

$$
a_\text{adv}^\text{int}(u, v) = -\sum_K \int_K (\mathbf{w} u) \cdot \nabla v\, \mathrm{d}x
+ \sum_{F \in \mathcal{F}_I} \int_F (\mathbf{w} \cdot \mathbf{n})\, u_\text{up}\, \llbracket v \rrbracket\, \mathrm{d}s
$$

Boundary contributions depend on the flow direction and are treated in the next section.

## Boundary Conditions

Each boundary condition modifies the face integrals that arise during integration by parts on $\mathcal{F}_B$.

### Diffusion

**Dirichlet** ($u = g_D$ on $\Gamma_D$):

Imposed *weakly* by treating $\Gamma_D$ faces as one-sided SIPG faces, mirroring the interior treatment:

$$
\begin{align}
a_\text{diff}^D(u, v) &= - \int_{\Gamma_D} D(\nabla u \cdot \mathbf{n})\, v\, \mathrm{d}s
  - \int_{\Gamma_D} D(\nabla v \cdot \mathbf{n})\, u\, \mathrm{d}s
  + \int_{\Gamma_D} \frac{\alpha D}{h}\, u\, v\, \mathrm{d}s \\
\ell_\text{diff}^D(v) &= - \int_{\Gamma_D} D(\nabla v \cdot \mathbf{n})\, g_D\, \mathrm{d}s
  + \int_{\Gamma_D} \frac{\alpha D}{h}\, g_D\, v\, \mathrm{d}s
\end{align}
$$

**Neumann** ($D\nabla u \cdot \mathbf{n} = g_N$ on $\Gamma_N$):

The boundary flux integral is replaced directly by the prescribed value. No penalty is needed:

$$
-\int_{\Gamma_N} D(\nabla u \cdot \mathbf{n})\, v\, \mathrm{d}s \longrightarrow -\int_{\Gamma_N} g_N\, v\, \mathrm{d}s
$$

The homogeneous case $g_N = 0$ is the natural condition and requires no modification.

**Robin** ($D\nabla u \cdot \mathbf{n} = g_R - \beta u$ on $\Gamma_R$, $\beta \geq 0$):

Substituting the Robin condition into the boundary face integral:

$$
-\int_{\Gamma_R} D(\nabla u \cdot \mathbf{n})\, v\, \mathrm{d}s = \int_{\Gamma_R} \beta\, u\, v\, \mathrm{d}s - \int_{\Gamma_R} g_R\, v\, \mathrm{d}s
$$

The $\beta u v$ term enters the bilinear form; $g_R v$ moves to the right-hand side.

### Advection

**Inflow** ($\mathbf{w} \cdot \mathbf{n} < 0$ on $\Gamma_\text{in}$, $u = g_\text{in}$):

The upwind value is the prescribed inflow data; information enters from outside the domain:

$$
\int_{\Gamma_\text{in}} (\mathbf{w} \cdot \mathbf{n})\, g_\text{in}\, v\, \mathrm{d}s \quad \text{(RHS only)}
$$

Since $\mathbf{w} \cdot \mathbf{n} < 0$ on $\Gamma_\text{in}$, this contributes positively for positive inflow data.

**Outflow** ($\mathbf{w} \cdot \mathbf{n} > 0$ on $\Gamma_\text{out}$):

The upwind value is the interior trace $u$; no data needs to be prescribed:

$$
\int_{\Gamma_\text{out}} (\mathbf{w} \cdot \mathbf{n})\, u\, v\, \mathrm{d}s \quad \text{(bilinear form)}
$$

This is the natural outflow condition for DG upwinding.

**Wall / Symmetry** ($\mathbf{w} \cdot \mathbf{n} = 0$ on $\Gamma_w$):

No advective flux crosses the boundary. The face integral vanishes, so no term is added for either walls or symmetry planes.

## The Complete Weak Formulation

Find $u \in V_h$ such that $a(u, v) = \ell(v)$ for all $v \in V_h$, where:

**Bilinear form**:

$$
\begin{split}
a(u, v)
&= \sum_K \int_K D\, \nabla u \cdot \nabla v\, \mathrm{d}x \\
&\quad - \sum_{F \in \mathcal{F}_I \cup \Gamma_D} \int_F \langle D\nabla u \rangle \cdot \mathbf{n}\, \llbracket v \rrbracket\, \mathrm{d}s
   - \sum_{F \in \mathcal{F}_I \cup \Gamma_D} \int_F \langle D\nabla v \rangle \cdot \mathbf{n}\, \llbracket u \rrbracket\, \mathrm{d}s
   + \sum_{F \in \mathcal{F}_I \cup \Gamma_D} \int_F \frac{\alpha D}{h}\, \llbracket u \rrbracket \llbracket v \rrbracket\, \mathrm{d}s \\
&\quad + \int_{\Gamma_R} \beta\, u\, v\, \mathrm{d}s \\
&\quad - \sum_K \int_K (\mathbf{w} u) \cdot \nabla v\, \mathrm{d}x
   + \sum_{F \in \mathcal{F}_I} \int_F (\mathbf{w} \cdot \mathbf{n})\, u_\text{up}\, \llbracket v \rrbracket\, \mathrm{d}s
   + \int_{\Gamma_\text{out}} (\mathbf{w} \cdot \mathbf{n})\, u\, v\, \mathrm{d}s
\end{split}
$$

**Linear form**:

$$
\begin{split}
\ell(v)
&= \int_\Omega f\, v\, \mathrm{d}x \\
&\quad - \int_{\Gamma_N} g_N\, v\, \mathrm{d}s
   - \int_{\Gamma_R} g_R\, v\, \mathrm{d}s \\
&\quad - \int_{\Gamma_D} D(\nabla v \cdot \mathbf{n})\, g_D\, \mathrm{d}s
   + \int_{\Gamma_D} \frac{\alpha D}{h}\, g_D\, v\, \mathrm{d}s \\
&\quad - \int_{\Gamma_\text{in}} (\mathbf{w} \cdot \mathbf{n})\, g_\text{in}\, v\, \mathrm{d}s
\end{split}
$$

## Coupling with an External Velocity Field

The formulation above assumes a velocity field that is either prescribed analytically or represented continuously, divergence-free in the discrete sense, and exactly tangential on walls. When $\mathbf{w}$ is instead supplied by an external solver such as OpenFOAM, each of these assumptions relaxes, and the face integrals need care. The issues divide into three groups.

### Single-valued face flux

OpenFOAM stores its velocity as cell-centred values together with a set of face-normal mass fluxes $\phi_F = \int_F (\rho \mathbf{w}) \cdot \mathbf{n}\,\mathrm{d}s$ that are, by construction, single-valued per face and conservative at the discrete level. Interpolating the cell-centred velocity into a DG space on the transport mesh and reconstructing $\mathbf{w}\cdot\mathbf{n}$ from it loses both properties: the two sides of an interior face disagree, and $\llbracket \mathbf{w}\cdot\mathbf{n} \rrbracket \ne 0$.

The compact identity used for the interior advection integral then no longer collapses to an upwind flux. Two equivalent fixes are available:

1. Replace the compact form by an explicit upwind flux that does not rely on single-valuedness of $\mathbf{w}\cdot\mathbf{n}$:
   $$
   \hat{f}_\text{adv} = \tfrac{1}{2} \bigl( \mathbf{w}_0\cdot\mathbf{n}\, u_0 + \mathbf{w}_1\cdot\mathbf{n}\, u_1 \bigr)
   + \tfrac{1}{2} \bigl| \langle \mathbf{w}\rangle\cdot\mathbf{n} \bigr|\, \llbracket u \rrbracket.
   $$
2. Transfer the face-flux field $\phi_F$ directly from OpenFOAM and use it as the advective transport coefficient on each face, replacing $\mathbf{w}\cdot\mathbf{n}$ by $\phi_F / (\rho\, |F|)$ (or by $\phi_F / |F|$ if a volumetric flux is used). The upwind branch is then selected on the sign of $\phi_F$. This is the option that preserves OpenFOAM's flux conservation and is generally preferred.

In either case the boundary face integrals $\int_{\Gamma_\text{in/out}}$ must use the same face flux, not the reconstructed $\mathbf{w}\cdot\mathbf{n}$, so that the inflow and outflow splits remain consistent.

### Discrete divergence

The physical field is divergence-free, but the interpolant $\mathbf{w}_h$ on the transport mesh is not, in general. The two forms

$$
\nabla\cdot(\mathbf{w} u) \quad\text{and}\quad \mathbf{w}\cdot\nabla u
$$

differ by $(\nabla\cdot\mathbf{w})\, u$, which is zero at the continuum level but acts as a spurious source when $\nabla\cdot\mathbf{w}_h \ne 0$. Two options are available:

- Use the skew-symmetric or non-conservative form $\mathbf{w}\cdot\nabla u$ when deriving the weak form. This trades local conservation for reduced sensitivity to interpolation error in $\mathbf{w}$ and is the pragmatic choice when the velocity is known only approximately.
- Project $\mathbf{w}$ onto a divergence-conforming space such as Raviart-Thomas or BDM, which preserves $\nabla\cdot\mathbf{w}_h = 0$ at the discrete level. This keeps the conservative form valid but adds a projection step and couples the transport mesh more tightly to the flow discretisation.

When the face-flux form (option 2 above) is used, conservation is enforced through $\phi_F$ directly and this divergence mismatch does not arise in the same way.

### Wall boundary flux

On no-slip or symmetry walls the physical condition is $\mathbf{w}\cdot\mathbf{n} = 0$, and the corresponding face integral vanishes. An interpolated velocity field rarely satisfies this exactly: a small residual normal component acts as a spurious inflow or outflow, with the sign and magnitude depending on the interpolation. The safest treatment is to enforce $\mathbf{w}\cdot\mathbf{n} = 0$ explicitly on wall facets, either by masking the integrand or by using the face-flux field $\phi_F$, which is identically zero on OpenFOAM wall patches and therefore introduces no leak.

### Mesh transfer

All of the above assumes that values defined on the OpenFOAM mesh have been transferred to the transport mesh. Mesh-to-mesh interpolation is itself a source of error, particularly at boundaries and in regions with large velocity gradients, and its treatment is outside the scope of this note. For the purposes of this formulation, we assume that either a face-flux field $\phi_F$ or a cell-centred velocity has been made available on the transport mesh by a method appropriate to the coupling.

## Further Reading

- Arnold et al. (2002), *Unified analysis of discontinuous Galerkin methods for elliptic problems*, SIAM J. Numer. Anal.
- Di Pietro & Ern (2012), *Mathematical Aspects of Discontinuous Galerkin Methods*, Springer
- [DOLFINx documentation](https://docs.fenicsproject.org/dolfinx/main/python/)
