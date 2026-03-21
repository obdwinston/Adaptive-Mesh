import csv
import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

# === CONFIGURATION ===

cx, cy = 5.0, 10.0
R = 0.5
U_inf = 1.0
Re = 40.0
n_sample = 37
probe_offset = 0.0
output_dir = "output"

# === DATA LOADING ===

vts_files = sorted(glob.glob(os.path.join(output_dir, "step_*_patch_*.vts")))
last_step = vts_files[-1].split("step_")[1].split("_patch_")[0]
last_patches = sorted(
    glob.glob(os.path.join(output_dir, f"step_{last_step}_patch_*.vts"))
)

best_mesh = None
best_dx = float("inf")
for f in last_patches:
    mesh = pv.read(f)
    xs = np.unique(mesh.points[:, 0])
    dx = xs[1] - xs[0] if len(xs) > 1 else float("inf")
    if dx < best_dx:
        best_dx = dx
        best_mesh = mesh

dx = best_dx
pts = best_mesh.points

# === PROBE POINTS ===

theta_deg = np.linspace(0, 180, n_sample)
theta_rad = np.radians(theta_deg)

probe_pts = np.zeros((n_sample, 3))
R_probe = R + probe_offset * dx
probe_pts[:, 0] = cx + R_probe * np.cos(theta_rad)
probe_pts[:, 1] = cy + R_probe * np.sin(theta_rad)

# === INTERPOLATION VISUALISATION ===

fig, ax = plt.subplots(figsize=(8, 8))

margin = 3 * dx
near = (np.abs(pts[:, 0] - cx) < R + margin) & (np.abs(pts[:, 1] - cy) < R + margin)
near_pts = pts[near]

xs_near = np.unique(near_pts[:, 0])
ys_near = np.unique(near_pts[:, 1])
for x in xs_near:
    ax.axvline(x - dx / 2, color="lightgray", linewidth=0.5, zorder=1)
ax.axvline(xs_near[-1] + dx / 2, color="lightgray", linewidth=0.5, zorder=1)
for y in ys_near:
    ax.axhline(y - dx / 2, color="lightgray", linewidth=0.5, zorder=1)
ax.axhline(ys_near[-1] + dx / 2, color="lightgray", linewidth=0.5, zorder=1)

ax.scatter(
    near_pts[:, 0], near_pts[:, 1], s=8, c="steelblue", zorder=3, label="Cell centres"
)
ax.add_patch(plt.Circle((cx, cy), R, fill=True, color="lightgray", alpha=0.3, zorder=2))
ax.add_patch(
    plt.Circle(
        (cx, cy),
        R,
        fill=False,
        color="black",
        linewidth=2,
        zorder=4,
        label="Cylinder surface",
    )
)
ax.scatter(
    probe_pts[:, 0],
    probe_pts[:, 1],
    s=40,
    c="red",
    marker="x",
    zorder=5,
    linewidths=1.5,
    label="Probe points",
)

px, py = probe_pts[10, 0], probe_pts[10, 1]
i_left = max(0, min(np.searchsorted(xs_near, px) - 1, len(xs_near) - 2))
j_low = max(0, min(np.searchsorted(ys_near, py) - 1, len(ys_near) - 2))

corners_x = [xs_near[i_left], xs_near[i_left + 1], xs_near[i_left + 1], xs_near[i_left]]
corners_y = [ys_near[j_low], ys_near[j_low], ys_near[j_low + 1], ys_near[j_low + 1]]

ax.scatter(
    corners_x,
    corners_y,
    s=80,
    c="orange",
    edgecolors="darkorange",
    linewidths=1.5,
    zorder=6,
    label="Interpolation stencil",
)
for cxi, cyi in zip(corners_x, corners_y):
    ax.plot(
        [px, cxi], [py, cyi], "--", color="orange", linewidth=1, alpha=0.7, zorder=4
    )
ax.scatter([px], [py], s=100, c="red", marker="x", zorder=7, linewidths=2.5)
ax.add_patch(
    plt.Polygon(
        [
            (xs_near[i_left] - dx / 2, ys_near[j_low] - dx / 2),
            (xs_near[i_left + 1] + dx / 2, ys_near[j_low] - dx / 2),
            (xs_near[i_left + 1] + dx / 2, ys_near[j_low + 1] + dx / 2),
            (xs_near[i_left] - dx / 2, ys_near[j_low + 1] + dx / 2),
        ],
        fill=True,
        facecolor="orange",
        alpha=0.1,
        edgecolor="orange",
        linewidth=2,
        zorder=2,
    )
)

ax.set_xlim(cx - R - margin, cx + R + margin)
ax.set_ylim(cy - R - margin, cy + R + margin)
ax.set_aspect("equal")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Grid Cells to Cylinder Surface (Bilinear Interpolation)")
ax.legend(loc="upper right", fontsize=9)
fig.tight_layout()
plt.savefig("verification_interpolation.png", dpi=150)
plt.show()

# === PRESSURE COMPARISON ===

result = pv.PolyData(probe_pts).sample(best_mesh)
cp_sim = result.point_data["pressure"] / (0.5 * U_inf**2)

theta_ref, cp_ref = [], []
with open("verification.csv", "r") as f:
    for row in csv.DictReader(f):
        theta_ref.append(float(row["theta"]))
        cp_ref.append(float(row["cp"]))

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(theta_deg, cp_sim, "o-", label="Current", markersize=3)
ax.plot(theta_ref, cp_ref, "s-", label="Reference", markersize=4, alpha=0.8)
ax.set_xlabel(r"$\theta$ (degrees from rear)")
ax.set_ylabel(r"$C_p$")
ax.set_title(f"Pressure Coefficient on Cylinder Surface (Re = {Re:.0f})")
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 180)
fig.tight_layout()
plt.savefig("verification_pressure.png", dpi=150)
plt.show()
