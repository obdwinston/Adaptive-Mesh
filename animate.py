"""Reads VTK output from the Fortran solver and produces a 2x2 animation.

Panels:
    top-left:     velocity magnitude with patch outlines and body
    top-right:    vorticity with patch outlines and body
    bottom-left:  pressure with patch outlines and body
    bottom-right: grid nodes coloured by refinement level
"""

import glob
import math
import os
import re
import sys

import matplotlib.animation as animation
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

try:
    import pyvista as pv
except ImportError:
    print("pyvista not found. Install with: pip install pyvista")
    sys.exit(1)


# === CONSTANTS ===

LEVEL_COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]
LEVEL_LABELS = ["level 0", "level 1", "level 2", "level 3", "level 4"]
RECT_COLORS = LEVEL_COLORS


# === FILE DISCOVERY ===


def discover_steps(output_dir="output"):
    """Find all unique step numbers in the output directory."""
    pattern = os.path.join(output_dir, "step_*_patch_*.vts")
    files = glob.glob(pattern)
    steps = sorted(set(int(re.search(r"step_(\d+)", f).group(1)) for f in files))
    return steps


# === DATA LOADING ===


def load_step(step, output_dir="output"):
    """Load all patch VTK files for a given step.

    Returns a list of dicts with keys: X, Y, speed, vorticity, pressure,
    level, nx, ny, x0, y0, x_end, y_end, dx, dy.
    """
    pattern = os.path.join(output_dir, f"step_{step:05d}_patch_*.vts")
    files = sorted(glob.glob(pattern))
    patches = []

    for f in files:
        mesh = pv.read(f)
        points = mesh.points
        nx_pts = mesh.dimensions[0]
        ny_pts = mesh.dimensions[1]

        X = points[:, 0].reshape(ny_pts, nx_pts)
        Y = points[:, 1].reshape(ny_pts, nx_pts)

        speed = mesh.point_data.get("speed", np.zeros(mesh.n_points))
        speed = speed.reshape(ny_pts, nx_pts)

        vorticity = mesh.point_data.get("vorticity", np.zeros(mesh.n_points))
        vorticity = vorticity.reshape(ny_pts, nx_pts)

        pressure = mesh.point_data.get("pressure", np.zeros(mesh.n_points))
        pressure = pressure.reshape(ny_pts, nx_pts)

        # infer spacing from point positions
        dx_val = float(X[0, 1] - X[0, 0]) if nx_pts > 1 else 1.0
        dy_val = float(Y[1, 0] - Y[0, 0]) if ny_pts > 1 else 1.0

        level = 0
        try:
            lv = mesh.field_data.get("level", None)
            if lv is not None:
                level = int(lv[0])
        except (KeyError, IndexError):
            pass

        patches.append(
            {
                "X": X,
                "Y": Y,
                "speed": speed,
                "vorticity": vorticity,
                "pressure": pressure,
                "level": int(level),
                "nx": nx_pts,
                "ny": ny_pts,
                "x0": float(X.min()),
                "y0": float(Y.min()),
                "x_end": float(X.max()),
                "y_end": float(Y.max()),
                "dx": float(dx_val),
                "dy": float(dy_val),
            }
        )

    # infer level from dx ratio when metadata is absent
    if patches:
        max_dx = max(p["dx"] for p in patches)
        for p in patches:
            if p["level"] == 0:
                ratio = max_dx / p["dx"]
                if ratio > 1.0 + 1e-6:
                    p["level"] = round(math.log(ratio) / math.log(2))

    return patches


# === ANIMATION ===


def animate_vtk(
    output_dir="output",
    save_as="animation.mp4",
    domain=(0.0, 0.0, 30.0, 20.0),
    bodies=None,
    body_file="body.dat",
    U_inf=1.0,
    fps=30,
):
    """Produce a 3x1 MP4 animation from VTK output files."""
    if bodies is None:
        bodies = [{"cx": 5.0, "cy": 10.0, "half_w": 0.5, "half_h": 0.5}]

    steps = discover_steps(output_dir)
    if not steps:
        print(f"No VTK files found in {output_dir}/")
        return

    print(f"Found {len(steps)} steps: {steps[0]} to {steps[-1]}")

    x_lo, y_lo, x_hi, y_hi = domain
    vmax = 1.5 * U_inf

    fig, axes = plt.subplots(3, 1, figsize=(10, 14))
    fig.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.03, hspace=0.15)

    div0 = make_axes_locatable(axes[0])
    cax_vel = div0.append_axes("right", size="2%", pad=0.1)
    div1 = make_axes_locatable(axes[1])
    cax_vor = div1.append_axes("right", size="2%", pad=0.1)
    div2 = make_axes_locatable(axes[2])
    cax_grid = div2.append_axes("right", size="2%", pad=0.1)
    cax_grid.axis("off")

    # load polygon from body file
    body_polygon = None
    if os.path.isfile(body_file):
        raw = np.loadtxt(body_file)
        if raw.ndim == 2 and raw.shape[1] >= 2:
            coords = raw[:, :2].copy()
            cx_raw = 0.5 * (coords[:, 0].min() + coords[:, 0].max())
            cy_raw = 0.5 * (coords[:, 1].min() + coords[:, 1].max())
            w_raw = coords[:, 0].max() - coords[:, 0].min()
            h_raw = coords[:, 1].max() - coords[:, 1].min()
            D_raw = max(w_raw, h_raw)
            body0 = bodies[0]
            scale = 2.0 * max(body0["half_w"], body0["half_h"]) / D_raw
            coords[:, 0] = (coords[:, 0] - cx_raw) * scale + body0["cx"]
            coords[:, 1] = (coords[:, 1] - cy_raw) * scale + body0["cy"]
            body_polygon = coords

    def draw_bodies(ax):
        """Draw the body polygon or fallback rectangles on the given axis."""
        if body_polygon is not None:
            poly = plt.Polygon(body_polygon, fc="grey", ec="black", lw=1.5, zorder=10)
            ax.add_patch(poly)
        else:
            for body in bodies:
                hw = body["half_w"]
                hh = body["half_h"]
                ax.add_patch(
                    mpatches.Rectangle(
                        (body["cx"] - hw, body["cy"] - hh),
                        2 * hw,
                        2 * hh,
                        fc="grey",
                        ec="black",
                        lw=1.5,
                        zorder=10,
                    )
                )

    def draw_patch_outlines(ax, patches):
        """Draw refinement-level rectangles around non-coarse patches."""
        seen = set()
        for p in patches:
            if p["level"] == 0:
                continue
            color = RECT_COLORS[p["level"]]
            label = LEVEL_LABELS[p["level"]] if p["level"] not in seen else None
            seen.add(p["level"])
            ax.add_patch(
                mpatches.Rectangle(
                    (p["x0"], p["y0"]),
                    p["x_end"] - p["x0"],
                    p["y_end"] - p["y0"],
                    fill=False,
                    edgecolor=color,
                    linewidth=1.5,
                    label=label,
                )
            )

    # create colorbars once
    wlim = max(vmax * 2, 1e-10)
    sm_vel = plt.cm.ScalarMappable(cmap="RdYlBu_r", norm=plt.Normalize(0, vmax))
    sm_vor = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-wlim, wlim))
    fig.colorbar(sm_vel, cax=cax_vel, label="Velocity")
    fig.colorbar(sm_vor, cax=cax_vor, label="Vorticity")

    def update(frame_idx):
        """Render one animation frame."""
        step = steps[frame_idx]
        patches = load_step(step, output_dir)

        for ax in axes:
            ax.cla()

        # top: velocity magnitude
        ax = axes[0]
        kw = dict(shading="auto", cmap="RdYlBu_r", vmin=0, vmax=vmax)
        for p in patches:
            ax.pcolormesh(p["X"], p["Y"], p["speed"], **kw)
        draw_bodies(ax)
        draw_patch_outlines(ax, patches)
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)
        ax.set_aspect("equal")
        if len(patches) > 1:
            ax.legend(loc="upper right", fontsize=6)

        # middle: vorticity
        ax = axes[1]
        for p in patches:
            ax.pcolormesh(
                p["X"],
                p["Y"],
                p["vorticity"],
                shading="auto",
                cmap="RdBu_r",
                vmin=-wlim,
                vmax=wlim,
            )
        draw_bodies(ax)
        draw_patch_outlines(ax, patches)
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)
        ax.set_aspect("equal")

        # bottom: grid nodes by level
        ax = axes[2]
        seen = set()
        for p in patches:
            color = LEVEL_COLORS[p["level"]]
            if p["level"] not in seen:
                label = (
                    f"{LEVEL_LABELS[p['level']]} (dx={p['dx']:.4f}, dy={p['dy']:.4f})"
                )
                seen.add(p["level"])
            else:
                label = None
            sz = max(0.3, 3.0 / (2 ** p["level"]))
            ax.scatter(
                p["X"].ravel(),
                p["Y"].ravel(),
                s=sz,
                color=color,
                label=label,
            )
        draw_bodies(ax)
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)
        ax.set_aspect("equal")
        ax.legend(loc="upper right", fontsize=5, markerscale=4)

    print(f"Rendering {len(steps)} frames...")
    anim = animation.FuncAnimation(
        fig,
        update,
        frames=len(steps),
        interval=200,
        repeat=False,
    )
    anim.save(save_as, writer="ffmpeg", fps=fps)
    print(f"Saved {save_as}")


if __name__ == "__main__":
    animate_vtk()
