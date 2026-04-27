"""
Illustrates why upwinding matters in advection-dominated problems.

Top row: cell-by-cell schematics showing how each scheme computes the face flux.
Bottom row: resulting 1D solutions compared against the exact answer.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec

# ── Parameters ────────────────────────────────────────────────────────────────
w = 1.0  # advective velocity (m/s)
D = 0.015  # diffusion coefficient (m²/s)
n_cells = 8  # number of cells
h = 1.0 / n_cells
Pe_cell = w * h / (2 * D)

# ── Node positions ────────────────────────────────────────────────────────────
nodes = np.linspace(0, 1, n_cells + 1)
cell_centres = 0.5 * (nodes[:-1] + nodes[1:])

# ── Exact solution ────────────────────────────────────────────────────────────
x_ex = np.linspace(0, 1, 500)
u_ex = (np.exp(w / D * x_ex) - 1) / (np.exp(w / D) - 1)

# ── CG solution (central differencing = standard linear CG) ──────────────────
# Equation at interior node i:
#   -(D/h² + w/2h) u_{i-1} + (2D/h²) u_i + (-D/h² + w/2h) u_{i+1} = 0
sz = n_cells - 1
A_cg = np.zeros((sz, sz))
for k in range(sz):
    A_cg[k, k] = 2 * D / h**2
    if k > 0:
        A_cg[k, k - 1] = -(D / h**2 + w / (2 * h))
    if k < sz - 1:
        A_cg[k, k + 1] = -D / h**2 + w / (2 * h)
b_cg = np.zeros(sz)
b_cg[-1] = D / h**2 - w / (2 * h)  # contribution from u_n = 1

u_cg = np.zeros(n_cells + 1)
u_cg[-1] = 1.0
u_cg[1:n_cells] = np.linalg.solve(A_cg, b_cg)

# ── Upwind solution (first-order upwind advection + central diffusion) ────────
# Equation at interior node i:
#   -(w/h + D/h²) u_{i-1} + (w/h + 2D/h²) u_i + (-D/h²) u_{i+1} = 0
A_up = np.zeros((sz, sz))
for k in range(sz):
    A_up[k, k] = w / h + 2 * D / h**2
    if k > 0:
        A_up[k, k - 1] = -(w / h + D / h**2)
    if k < sz - 1:
        A_up[k, k + 1] = -D / h**2
b_up = np.zeros(sz)
b_up[-1] = D / h**2  # contribution from u_n = 1

u_up = np.zeros(n_cells + 1)
u_up[-1] = 1.0
u_up[1:n_cells] = np.linalg.solve(A_up, b_up)

# Cell averages used to colour the schematic cells
u_cg_cells = 0.5 * (u_cg[:-1] + u_cg[1:])
u_up_cells = 0.5 * (u_up[:-1] + u_up[1:])

# ── Colours ───────────────────────────────────────────────────────────────────
COL_CG = "#4dac26"
COL_DG = "#2166ac"
COL_FLOW = "#d62728"


# ── Schematic drawing ─────────────────────────────────────────────────────────
def draw_schematic(ax, u_cells, upwind, title, arrow_col):
    """
    Draw n_show cells as coloured rectangles.
    At each interior face show which cell(s) contribute to the face flux:
      upwind=False  ->  arrows from both neighbours (CG / symmetric)
      upwind=True   ->  arrow from upstream cell only (DG / upwind)
    """
    n_show = 5
    cw, ch = 1.2, 0.55  # cell width / height

    cmap = plt.cm.Blues
    u_clamp = np.clip(u_cells[:n_show], 0.0, 1.0)

    # Draw cells
    for i in range(n_show):
        x0 = i * cw
        rect = mpatches.FancyBboxPatch(
            [x0, 0],
            cw,
            ch,
            boxstyle="square,pad=0",
            linewidth=2,
            edgecolor="#333333",
            facecolor=cmap(0.15 + 0.8 * u_clamp[i]),
        )
        ax.add_patch(rect)
        txt_col = "white" if u_clamp[i] > 0.55 else "black"
        ax.text(
            x0 + cw / 2,
            ch / 2,
            f"{u_cells[i]:.2f}",
            ha="center",
            va="center",
            fontsize=11,
            color=txt_col,
            fontweight="bold",
        )

    # Flow direction arrow
    ax.annotate(
        "",
        xy=(n_show * cw - 0.15, ch + 0.18),
        xytext=(0.15, ch + 0.18),
        arrowprops=dict(arrowstyle="->", color=COL_FLOW, lw=2.2),
    )
    ax.text(
        n_show * cw / 2,
        ch + 0.31,
        "flow direction",
        ha="center",
        fontsize=9,
        color=COL_FLOW,
    )

    # Face flux arrows
    for i in range(1, n_show):
        xf = i * cw
        yf = ch / 2
        if upwind:
            # Single arrow: information comes from upstream (left) cell only
            ax.annotate(
                "",
                xy=(xf + 0.07, yf),
                xytext=(xf - 0.32, yf),
                arrowprops=dict(arrowstyle="->", color=arrow_col, lw=1.8),
            )
        else:
            # Two arrows: both neighbours contribute (symmetric)
            ax.annotate(
                "",
                xy=(xf, yf + 0.07),
                xytext=(xf - 0.32, yf + 0.07),
                arrowprops=dict(arrowstyle="->", color=arrow_col, lw=1.8),
            )
            ax.annotate(
                "",
                xy=(xf, yf - 0.07),
                xytext=(xf + 0.32, yf - 0.07),
                arrowprops=dict(arrowstyle="->", color=arrow_col, lw=1.8),
            )

    lbl = "upstream value only" if upwind else "average of both neighbours"
    ax.text(
        n_show * cw / 2,
        -0.17,
        f"Face flux uses: {lbl}",
        ha="center",
        fontsize=9,
        color=arrow_col,
        style="italic",
    )

    ax.set_xlim(-0.25, n_show * cw + 0.25)
    ax.set_ylim(-0.22, ch + 0.42)
    ax.axis("off")
    ax.set_title(title, fontsize=11, fontweight="bold", pad=4)


# ── Build figure ──────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(13, 8))
gs = gridspec.GridSpec(2, 2, height_ratios=[0.55, 1.8], hspace=0.30, wspace=0.35)
fig.subplots_adjust(top=0.93, bottom=0.09, left=0.06, right=0.97)

ax_cg_s = fig.add_subplot(gs[0, 0])
ax_dg_s = fig.add_subplot(gs[0, 1])
ax_cg_p = fig.add_subplot(gs[1, 0])
ax_dg_p = fig.add_subplot(gs[1, 1])

draw_schematic(
    ax_cg_s, u_cg_cells, upwind=False, title="CG: symmetric face flux", arrow_col=COL_CG
)

draw_schematic(
    ax_dg_s, u_up_cells, upwind=True, title="DG: upwind face flux", arrow_col=COL_DG
)

# Solution plots
for ax, u_num, col, title in [
    (ax_cg_p, u_cg, COL_CG, "CG solution"),
    (ax_dg_p, u_up, COL_DG, "DG / upwind solution"),
]:
    ax.plot(x_ex, u_ex, "k--", lw=1.5, label="Exact", zorder=3)
    ax.plot(
        nodes, u_num, color=col, lw=1.8, marker="o", ms=5, label="Numerical", zorder=4
    )
    ax.axhline(0, color="grey", lw=0.7, ls=":")
    ax.axhline(1, color="grey", lw=0.7, ls=":")
    ax.set_xlabel("Position $x$", fontsize=11)
    ax.set_ylabel("Concentration $u$", fontsize=11)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, loc="upper left")
    ax.set_xlim(0, 1)
    ax.grid(True, alpha=0.25)
    ax.text(
        0.97,
        0.08,
        rf"Pe$_\mathrm{{cell}}$ = {Pe_cell:.1f}",
        transform=ax.transAxes,
        fontsize=9,
        ha="right",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.6),
    )

ax_cg_p.set_ylim(-0.5, 1.5)

fig.suptitle(
    "Upwinding: how the choice of face flux affects stability at high Péclet number",
    fontsize=12,
    fontweight="bold",
)

plt.savefig("images/cg_vs_dg_upwinding.svg", bbox_inches="tight")
plt.savefig("images/cg_vs_dg_upwinding.png", bbox_inches="tight", dpi=150)
plt.show()

print(f"Cell Peclet number : {Pe_cell:.2f}")
print(f"CG solution range  : [{u_cg.min():.3f}, {u_cg.max():.3f}]")
print(f"Upwind range       : [{u_up.min():.3f}, {u_up.max():.3f}]")
