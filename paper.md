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

Advection-diffusion equations describe how a quantity is transported through a medium by two competing mechanisms: advection (bulk transport carried along by a flow) and diffusion (spreading down a concentration gradient). In this work the quantity of interest is a species concentration in a fluid, specifically hydrogen, but the same equation governs heat transport, pollutant dispersion, and many other physical problems.

We consider the steady-state form for a scalar quantity $u$ (m$^{-3}$):

$$
\nabla \cdot (\mathbf{w} u) - \nabla \cdot (D \nabla u) = f \quad \text{in } \Omega
$$

where $\mathbf{w}$ is a prescribed advective velocity field (m$\,$s$^{-1}$), $D$ is the diffusion coefficient (m$^{2}\,$s$^{-1}$), $f$ is a volumetric source term (m$^{-3}\,$s$^{-1}$), and $\Omega \subset \mathbb{R}^d$ is the computational domain.

The relative importance of the two transport mechanisms is captured by the cell Péclet number:

$$
\text{Pe} = \frac{|\mathbf{w}| h}{2D}
$$

where $h$ is the local mesh element size. The Péclet number is a dimensionless ratio: large $|\mathbf{w}|$ or small $D$ pushes it up, and small $|\mathbf{w}|$ or large $D$ pushes it down.

- When $\text{Pe} \ll 1$, diffusion dominates. The solution is smooth and the problem behaves like a standard Poisson (diffusion) equation. This type of problem is called *elliptic*, meaning the solution at any point is influenced smoothly by conditions everywhere in the domain, rather than being driven by information propagating in a particular direction.
- When $\text{Pe} \gg 1$, advection dominates. The solution can develop *sharp layers*, which are thin regions of rapid concentration change near inflow boundaries or obstacles. These are physically real but numerically difficult to resolve accurately.

## Numerical Challenges and the Case for DG

Standard continuous Galerkin (CG) finite element methods work well when diffusion dominates, but as the Péclet number grows the solution develops oscillations that are entirely numerical in origin and have no physical meaning. To understand why, it helps to think about how information travels in an advection-dominated problem. If a fluid is moving to the right, the concentration at a point is determined by what is upstream, not downstream. A numerical scheme that respects this directionality is said to have *upwinding*: it biases the treatment of the advection term towards information coming from the upstream direction. Standard CG treats information from both sides of an element face symmetrically, which is appropriate for diffusion but wrong for advection, and this mismatch is what produces the spurious oscillations.

Several stabilisation strategies have been developed to recover stability within the CG framework:

- **Isotropic artificial diffusion** adds an extra diffusion term $D_\text{art} = \delta h |\mathbf{w}|$ uniformly in all directions, where $\delta$ is a dimensionless tuning parameter and $h$ is the local mesh size. It is simple to implement but *inconsistent*: the added term does not vanish even when the exact solution is substituted, meaning it permanently modifies the problem being solved. This introduces cross-stream smearing and reduces accuracy, particularly near sharp layers [(COMSOL, 2020)](https://www.comsol.com/blogs/understanding-stabilization-methods).
- **Streamline-upwind Petrov–Galerkin (SUPG)** and **Galerkin least-squares (GLS)** are smarter approaches that add numerical diffusion only along the flow direction, not across it. They are *consistent* stabilisations, meaning the extra terms vanish exactly when the true solution is substituted into the equations, so the formal accuracy of the method is preserved. These are generally preferred over isotropic diffusion.

These methods can be effective, but all require tuning parameters, add complexity to the formulation, and are fundamentally workarounds for a framework not designed with advection-dominated transport in mind.

DG methods address the root cause rather than patching the symptom. In DG the approximate solution is allowed to be discontinuous across element boundaries, and information is passed between elements through *numerical flux functions* defined on each shared face. The choice of flux naturally encodes upwinding for the advection term, without any add-on stabilisation. This gives DG several structural advantages:

- **Natural upwinding**: upwind information transfer is built directly into the flux formulation, so no stabilisation parameter needs to be tuned for the advection term.
- **Local conservation**: the flux balance is satisfied element by element, which is important for transport quantities.
- **Flexible boundary conditions**: all conditions are imposed weakly through face integrals, giving a consistent treatment across different condition types.
- **High-order accuracy**: high-degree polynomial spaces can be used on unstructured meshes without additional complications.

The trade-off is a larger global system (DG has more degrees of freedom than CG for the same mesh) and the need to choose a penalty parameter for the diffusion term carefully. Both are manageable in practice.

## Derivation of the DG Weak Form

### Notation

The domain $\Omega$ is divided into a *mesh* $\mathcal{T}_h = \{K\}$: a collection of non-overlapping elements $K$, triangles in 2D or tetrahedra in 3D, that together cover $\Omega$. The subscript $h$ refers to the typical element size.

Within this mesh, every shared boundary between two adjacent elements is an *interior face*. The set of all interior faces is $\mathcal{F}_I$, and the set of faces on the outer boundary of $\Omega$ is $\mathcal{F}_B$.

In DG, the approximate solution is allowed to be a different polynomial on each element and may be discontinuous across element boundaries. The *function space* that encodes this is:

$$
V_h = \{v \in L^2(\Omega) : v|_K \in \mathbb{P}_p(K),\ \forall K \in \mathcal{T}_h\}
$$

This says: $V_h$ contains all functions that, when restricted to any element $K$, are a polynomial of degree at most $p$, with no continuity requirement between elements.

Because solutions in $V_h$ can be discontinuous, each interior face $F \in \mathcal{F}_I$ has *two* values, one from each neighbouring element $K_0$ and $K_1$. Let $\mathbf{n} = \mathbf{n}_0$ be the outward unit normal from $K_0$ (pointing into $K_1$). For a quantity $q$ with values $q_0$ and $q_1$ on each side, we define the **average** and **jump** operators:

$$
\langle q \rangle := \frac{q_0 + q_1}{2} \quad \text{(average)}, \qquad \llbracket q \rrbracket := q_0 - q_1 \quad \text{(jump)}
$$

The average is the mean value across the face; the jump measures the size of the discontinuity, which is zero if the solution is continuous at that face. For a vector $\boldsymbol{\tau}$:

$$
\langle \boldsymbol{\tau} \rangle := \frac{\boldsymbol{\tau}_0 + \boldsymbol{\tau}_1}{2}, \qquad \llbracket \boldsymbol{\tau} \rrbracket := \boldsymbol{\tau}_0 \cdot \mathbf{n}_0 + \boldsymbol{\tau}_1 \cdot \mathbf{n}_1
$$

where $\llbracket \boldsymbol{\tau} \rrbracket = (\boldsymbol{\tau}_0 - \boldsymbol{\tau}_1) \cdot \mathbf{n}$ is the scalar normal jump. On a boundary face there is only one element, so $\langle q \rangle = q$ and $\llbracket q \rrbracket = q$.

### From Strong Form to Weak Form

The *strong form* of the PDE is an equation that must hold pointwise everywhere in $\Omega$, requiring $u$ to be smooth enough to be differentiated twice. On a mesh where the solution can be discontinuous, this is too strict.

The *weak form* relaxes this by multiplying the PDE by a *test function* $v$ and integrating over the domain. Rather than demanding the equation hold at every point, we demand it holds in an average sense against every possible $v$. The test function is chosen from the same space $V_h$ as the solution, which is called the *Galerkin* approach.

We derive the weak form element by element. Multiplying by $v$ and integrating over a single element $K$ gives:

$$
\int_K \nabla \cdot (\mathbf{w} u)\, v\, \mathrm{d}x - \int_K \nabla \cdot (D\nabla u)\, v\, \mathrm{d}x = \int_K f\, v\, \mathrm{d}x
$$

Applying the *divergence theorem* (the multi-dimensional analogue of integration by parts) moves one derivative off $u$ and onto $v$, reducing the smoothness requirement on $u$. It also produces boundary integrals on $\partial K$, which represent the flux of information through the faces of element $K$:

$$
\begin{align}
&\underbrace{-\int_K (\mathbf{w} u) \cdot \nabla v\, \mathrm{d}x}_{\text{advection, bulk}}
+ \underbrace{\int_{\partial K} (\mathbf{w} \cdot \mathbf{n}_K)\, u\, v\, \mathrm{d}s}_{\text{advection, face flux}} \\
&+ \underbrace{\int_K D\nabla u \cdot \nabla v\, \mathrm{d}x}_{\text{diffusion, bulk}}
- \underbrace{\int_{\partial K} D(\nabla u \cdot \mathbf{n}_K)\, v\, \mathrm{d}s}_{\text{diffusion, face flux}}
= \int_K f\, v\, \mathrm{d}x
\end{align}
$$

To form the global problem we sum over all elements $K \in \mathcal{T}_h$. The face integrals on $\partial K$ are shared between neighbouring elements and must be assembled carefully using the average and jump operators. The treatment of the advection and diffusion face terms differs, and is described in the two sections below.

### Diffusion Term: SIPG

When we sum the diffusion face integrals over all elements, each interior face $F$ is counted twice, once from each side. Using the average and jump operators this collapses to:

$$
\sum_K \int_{\partial K} D(\nabla u \cdot \mathbf{n}_K)\, v\, \mathrm{d}s
= \sum_{F \in \mathcal{F}_I} \int_F \left( \langle D\nabla u \rangle \cdot \mathbf{n}\, \llbracket v \rrbracket + \llbracket D\nabla u \rrbracket\, \langle v \rangle \right) \mathrm{d}s
+ \sum_{F \in \mathcal{F}_B} \int_F D(\nabla u \cdot \mathbf{n})\, v\, \mathrm{d}s
$$

For a smooth exact solution, the normal gradient is the same on both sides of every face, so $\llbracket D\nabla u \rrbracket = 0$ and the second interior term vanishes. We can drop it without changing which solutions satisfy the equations, a property called *consistency*. After dropping it, however, the resulting form is not symmetric (i.e. $a(u,v) \ne a(v,u)$) and is not yet stable. The Symmetric Interior Penalty Galerkin (SIPG) method adds two further face terms to fix both issues:

$$
\begin{align}
a_\text{diff}(u, v)
&= \underbrace{\sum_K \int_K D\, \nabla u \cdot \nabla v\, \mathrm{d}x}_{\text{bulk diffusion}} \\
&\quad - \underbrace{\sum_{F \in \mathcal{F}_I} \int_F \langle D\nabla u \rangle \cdot \mathbf{n}\, \llbracket v \rrbracket\, \mathrm{d}s}_{\text{consistency}} \\
&\quad - \underbrace{\sum_{F \in \mathcal{F}_I} \int_F \langle D\nabla v \rangle \cdot \mathbf{n}\, \llbracket u \rrbracket\, \mathrm{d}s}_{\text{symmetry}} \\
&\quad + \underbrace{\sum_{F \in \mathcal{F}_I} \int_F \frac{\alpha D}{h}\, \llbracket u \rrbracket \llbracket v \rrbracket\, \mathrm{d}s}_{\text{penalty}}
\end{align}
$$

The role of each term:

- **Bulk diffusion**: the standard $\int D\nabla u \cdot \nabla v$ integral, identical to standard Galerkin for pure diffusion.
- **Consistency**: the face flux term from integration by parts. It vanishes when the exact solution is substituted, so it does not change which solutions the equations accept.
- **Symmetry**: the consistency term with $u$ and $v$ swapped. It is also zero for the exact solution, but its presence makes $a(u,v) = a(v,u)$, which is required for a well-posed system.
- **Penalty**: penalises the size of the jump $\llbracket u \rrbracket$ across interior faces. This prevents the solution from becoming arbitrarily discontinuous and ensures the method is *coercive*, meaning the bilinear form is bounded away from zero, which guarantees stability. The penalty parameter is $\alpha = Cp^2$ (with $C \sim 10$ and $p$ the polynomial degree); it must be chosen large enough for stability, but not so large that it over-constrains the solution.

### Advection Term: Upwind Flux

When the advection face integrals are summed over all elements, each interior face $F$ produces a term involving $u$ evaluated from both sides. Since $u$ is discontinuous, we must choose which value to use, and this is where upwinding is applied.

The *upwind value* of $u$ on face $F$ is:

$$
u_\text{up} = \begin{cases} u_0 & \text{if } \mathbf{w} \cdot \mathbf{n} > 0 \quad (\text{flow from } K_0 \text{ to } K_1) \\ u_1 & \text{if } \mathbf{w} \cdot \mathbf{n} \leq 0 \quad (\text{flow from } K_1 \text{ to } K_0) \end{cases}
$$

Always take $u$ from the element the flow is coming *from*. This is the DG implementation of the upwinding principle, and it can be written compactly as:

$$
(\mathbf{w} \cdot \mathbf{n})\, u_\text{up} = \langle \mathbf{w} u \rangle \cdot \mathbf{n} + \tfrac{1}{2}|\mathbf{w} \cdot \mathbf{n}|\, \llbracket u \rrbracket
$$

The interior advection contribution is then:

$$
\begin{align}
a_\text{adv}(u, v)
&= \underbrace{-\sum_K \int_K (\mathbf{w} u) \cdot \nabla v\, \mathrm{d}x}_{\text{bulk advection}} \\
&\quad + \underbrace{\sum_{F \in \mathcal{F}_I} \int_F (\mathbf{w} \cdot \mathbf{n})\, u_\text{up}\, \llbracket v \rrbracket\, \mathrm{d}s}_{\text{interior upwind flux}}
\end{align}
$$

The **bulk advection** term moves information within each element; the **interior upwind flux** passes information between elements using the upstream value. Together they give DG its natural stability for advection without any additional tuning parameter.

Note: the compact identity for $u_\text{up}$ assumes $\mathbf{w}\cdot\mathbf{n}$ is single-valued on $F$. This holds when $\mathbf{w}$ is prescribed analytically or lives in a continuous FE space, but requires care when the velocity is interpolated from an external solver, as discussed in the section on coupling with an external velocity field.

Boundary contributions depend on the flow direction and are treated in the next section.

## Boundary Conditions

Each type of boundary condition enters the formulation by modifying the face integrals on $\mathcal{F}_B$ that appeared during integration by parts. On a boundary face there is only one element, so there is no "other side" to take values from. Instead, the boundary data closes the equations.

### Diffusion

**Dirichlet** ($u = g_D$ on $\Gamma_D$):

A Dirichlet condition prescribes the solution value directly. In a classical finite element method this would mean enforcing $u = g_D$ at shared nodes. In DG, where there are no shared nodes between elements, the condition is instead imposed *weakly* through the face integrals, by treating the boundary face exactly like an interior SIPG face but with $g_D$ playing the role of the solution on the outside:

$$
\begin{align}
a_\text{diff}^D(u, v)
&= \underbrace{- \int_{\Gamma_D} D(\nabla u \cdot \mathbf{n})\, v\, \mathrm{d}s}_{\text{flux consistency}} \\
&\quad - \underbrace{\int_{\Gamma_D} D(\nabla v \cdot \mathbf{n})\, u\, \mathrm{d}s}_{\text{symmetry}} \\
&\quad + \underbrace{\int_{\Gamma_D} \frac{\alpha D}{h}\, u\, v\, \mathrm{d}s}_{\text{penalty}} \\[6pt]
\ell_\text{diff}^D(v)
&= \underbrace{- \int_{\Gamma_D} D(\nabla v \cdot \mathbf{n})\, g_D\, \mathrm{d}s}_{\text{symmetry (data)}} \\
&\quad + \underbrace{\int_{\Gamma_D} \frac{\alpha D}{h}\, g_D\, v\, \mathrm{d}s}_{\text{penalty (data)}}
\end{align}
$$

The bilinear form $a_\text{diff}^D$ contains terms involving the unknown $u$ and goes on the left-hand side; the linear form $\ell_\text{diff}^D$ contains the known data $g_D$ and goes on the right-hand side.

**Neumann** ($D\nabla u \cdot \mathbf{n} = g_N$ on $\Gamma_N$):

A Neumann condition prescribes the flux through the boundary rather than the value. The diffusion face flux integral that arose from integration by parts is simply replaced by the prescribed value:

$$
-\int_{\Gamma_N} D(\nabla u \cdot \mathbf{n})\, v\, \mathrm{d}s \longrightarrow -\int_{\Gamma_N} g_N\, v\, \mathrm{d}s
$$

No penalty term is needed because we are specifying the gradient, not enforcing a value of $u$. The homogeneous case $g_N = 0$ (zero flux, such as an insulating or symmetry boundary) is the *natural* boundary condition for this formulation: it requires no additional terms at all, and is automatically satisfied if the boundary face integrals are simply omitted.

**Robin** ($D\nabla u \cdot \mathbf{n} = g_R - \beta u$ on $\Gamma_R$, $\beta \geq 0$):

A Robin condition is a linear combination of the flux and the solution value. It arises in heat transfer (convective cooling at a surface) and in mass transfer (surface recombination of species). Substituting the Robin condition into the flux integral splits it into one term depending on $u$ and one depending on the data $g_R$:

$$
-\int_{\Gamma_R} D(\nabla u \cdot \mathbf{n})\, v\, \mathrm{d}s = \underbrace{\int_{\Gamma_R} \beta\, u\, v\, \mathrm{d}s}_{\text{bilinear form}} - \underbrace{\int_{\Gamma_R} g_R\, v\, \mathrm{d}s}_{\text{linear form}}
$$

### Advection

The advection term also produces boundary face integrals, and their treatment depends on whether flow is entering or leaving the domain at each face.

**Inflow** ($\mathbf{w} \cdot \mathbf{n} < 0$ on $\Gamma_\text{in}$):

The flow is entering the domain, so the upwind value lies *outside* the domain. We prescribe it as $g_\text{in}$, which moves the term to the right-hand side:

$$
-\int_{\Gamma_\text{in}} (\mathbf{w} \cdot \mathbf{n})\, g_\text{in}\, v\, \mathrm{d}s \quad \text{(linear form)}
$$

The minus sign and the fact that $\mathbf{w} \cdot \mathbf{n} < 0$ on $\Gamma_\text{in}$ mean this contributes positively for positive inflow data, which is physically correct: inflow adds concentration to the domain.

**Outflow** ($\mathbf{w} \cdot \mathbf{n} > 0$ on $\Gamma_\text{out}$):

The flow is leaving the domain, so the upwind value is the interior trace $u$. No data needs to be prescribed; the term enters the bilinear form:

$$
\int_{\Gamma_\text{out}} (\mathbf{w} \cdot \mathbf{n})\, u\, v\, \mathrm{d}s \quad \text{(bilinear form)}
$$

**Wall / Symmetry** ($\mathbf{w} \cdot \mathbf{n} = 0$ on $\Gamma_w$):

No advective flux crosses the boundary. The face integral is identically zero, so no term needs to be added.

## The Complete Weak Formulation

Collecting all contributions from the interior and boundary face integrals, the complete problem is: find $u \in V_h$ such that $a(u, v) = \ell(v)$ for all $v \in V_h$, where the *bilinear form* $a$ collects all terms involving the unknown $u$ and the *linear form* $\ell$ collects all terms involving known data.

**Bilinear form**:

$$
\begin{align}
a(u, v)
&= \underbrace{\sum_K \int_K D\, \nabla u \cdot \nabla v\, \mathrm{d}x}_{\text{bulk diffusion}} \\
&\quad - \underbrace{\sum_{F \in \mathcal{F}_I \cup \Gamma_D} \int_F \langle D\nabla u \rangle \cdot \mathbf{n}\, \llbracket v \rrbracket\, \mathrm{d}s}_{\text{diffusion consistency}} \\
&\quad - \underbrace{\sum_{F \in \mathcal{F}_I \cup \Gamma_D} \int_F \langle D\nabla v \rangle \cdot \mathbf{n}\, \llbracket u \rrbracket\, \mathrm{d}s}_{\text{diffusion symmetry}} \\
&\quad + \underbrace{\sum_{F \in \mathcal{F}_I \cup \Gamma_D} \int_F \frac{\alpha D}{h}\, \llbracket u \rrbracket \llbracket v \rrbracket\, \mathrm{d}s}_{\text{diffusion penalty}} \\
&\quad + \underbrace{\int_{\Gamma_R} \beta\, u\, v\, \mathrm{d}s}_{\text{Robin}} \\
&\quad - \underbrace{\sum_K \int_K (\mathbf{w} u) \cdot \nabla v\, \mathrm{d}x}_{\text{bulk advection}} \\
&\quad + \underbrace{\sum_{F \in \mathcal{F}_I} \int_F (\mathbf{w} \cdot \mathbf{n})\, u_\text{up}\, \llbracket v \rrbracket\, \mathrm{d}s}_{\text{interior upwind flux}} \\
&\quad + \underbrace{\int_{\Gamma_\text{out}} (\mathbf{w} \cdot \mathbf{n})\, u\, v\, \mathrm{d}s}_{\text{outflow BC}}
\end{align}
$$

**Linear form**:

$$
\begin{align}
\ell(v)
&= \underbrace{\int_\Omega f\, v\, \mathrm{d}x}_{\text{source}} \\
&\quad - \underbrace{\int_{\Gamma_N} g_N\, v\, \mathrm{d}s}_{\text{Neumann BC}} \\
&\quad - \underbrace{\int_{\Gamma_R} g_R\, v\, \mathrm{d}s}_{\text{Robin BC}} \\
&\quad - \underbrace{\int_{\Gamma_D} D(\nabla v \cdot \mathbf{n})\, g_D\, \mathrm{d}s}_{\text{Dirichlet symmetry}} \\
&\quad + \underbrace{\int_{\Gamma_D} \frac{\alpha D}{h}\, g_D\, v\, \mathrm{d}s}_{\text{Dirichlet penalty}} \\
&\quad - \underbrace{\int_{\Gamma_\text{in}} (\mathbf{w} \cdot \mathbf{n})\, g_\text{in}\, v\, \mathrm{d}s}_{\text{inflow BC}}
\end{align}
$$

## Coupling with an External Velocity Field

The formulation above assumes the velocity field $\mathbf{w}$ is either prescribed analytically or lives in a continuous finite element space, is divergence-free (i.e. $\nabla\cdot\mathbf{w} = 0$) in the discrete sense, and has exactly zero normal component on walls. When $\mathbf{w}$ is supplied by an external solver such as OpenFOAM, each of these assumptions requires care. The issues divide into three groups.

### Single-valued face flux

In the interior advection term, we used the compact identity

$$
(\mathbf{w} \cdot \mathbf{n})\, u_\text{up} = \langle \mathbf{w} u \rangle \cdot \mathbf{n} + \tfrac{1}{2}|\mathbf{w} \cdot \mathbf{n}|\, \llbracket u \rrbracket
$$

which assumes $\mathbf{w}\cdot\mathbf{n}$ is the same on both sides of every interior face, i.e. $\llbracket \mathbf{w}\cdot\mathbf{n} \rrbracket = 0$. OpenFOAM satisfies this naturally: it stores face-normal mass fluxes $\phi_F = \int_F (\rho \mathbf{w}) \cdot \mathbf{n}\,\mathrm{d}s$ that are single-valued per face by construction. However, if the cell-centred velocity is interpolated into a DG space on the transport mesh and $\mathbf{w}\cdot\mathbf{n}$ is reconstructed from it, the two sides of each interior face will in general give different values, breaking the identity. Two fixes are available:

1. Replace the compact identity with an explicit upwind flux that does not assume single-valuedness:
   $$
   \hat{f}_\text{adv} = \tfrac{1}{2} \bigl( \mathbf{w}_0\cdot\mathbf{n}\, u_0 + \mathbf{w}_1\cdot\mathbf{n}\, u_1 \bigr)
   + \tfrac{1}{2} \bigl| \langle \mathbf{w}\rangle\cdot\mathbf{n} \bigr|\, \llbracket u \rrbracket
   $$
2. Transfer the face-flux field $\phi_F$ directly from OpenFOAM and use it as the advective transport coefficient on each face, replacing $\mathbf{w}\cdot\mathbf{n}$ by $\phi_F / (\rho\, |F|)$ (or by $\phi_F / |F|$ for a volumetric flux). The upwind branch is selected on the sign of $\phi_F$. This option preserves OpenFOAM's flux conservation exactly and is generally preferred.

In either case the boundary face integrals on $\Gamma_\text{in}$ and $\Gamma_\text{out}$ must use the same face flux consistently, not a reconstructed $\mathbf{w}\cdot\mathbf{n}$, so that the inflow and outflow classification remains consistent with the interior treatment.

### Discrete divergence

A *divergence-free* velocity field satisfies $\nabla\cdot\mathbf{w} = 0$, meaning there is no net creation or destruction of fluid volume at any point. The physical flow is divergence-free, but an interpolated velocity field $\mathbf{w}_h$ on the transport mesh generally is not. This matters because the two advection forms

$$
\nabla\cdot(\mathbf{w} u) \quad\text{and}\quad \mathbf{w}\cdot\nabla u
$$

differ by $(\nabla\cdot\mathbf{w})\, u$, which is zero in the continuous problem but acts as a spurious source when $\nabla\cdot\mathbf{w}_h \ne 0$. Two options are available:

- Use the non-conservative form $\mathbf{w}\cdot\nabla u$ when deriving the weak form. This reduces sensitivity to interpolation error in $\mathbf{w}$ at the cost of losing element-by-element conservation, and is the pragmatic choice when the velocity is known only approximately.
- Project $\mathbf{w}$ onto a divergence-conforming finite element space (such as Raviart-Thomas or BDM), which enforces $\nabla\cdot\mathbf{w}_h = 0$ exactly at the discrete level. This keeps the conservative form valid but adds a projection step and couples the transport mesh more tightly to the flow discretisation.

When the face-flux form (option 2 above) is used, conservation is enforced through $\phi_F$ directly and this divergence issue does not arise in the same way.

### Wall boundary flux

On no-slip or symmetry walls the physical condition is $\mathbf{w}\cdot\mathbf{n} = 0$, and the corresponding advection face integral vanishes. An interpolated velocity field rarely satisfies this exactly, and a small residual normal component acts as a spurious inflow or outflow. The safest fix is to enforce $\mathbf{w}\cdot\mathbf{n} = 0$ explicitly on wall faces, either by zeroing the integrand or by using the face-flux field $\phi_F$, which is identically zero on OpenFOAM wall patches and therefore introduces no spurious flux.

### Mesh transfer

All of the above assumes that velocity data defined on the OpenFOAM mesh has already been transferred to the transport mesh. Mesh-to-mesh interpolation is itself a source of error, particularly near boundaries and in regions with large velocity gradients, and its treatment is outside the scope of this document. We assume that either a face-flux field $\phi_F$ or a cell-centred velocity field has been made available on the transport mesh by an appropriate interpolation method.

## Further Reading
