#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import sys
import warnings
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")  # headless rendering – no display required
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.optimize import least_squares

try:
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
except ImportError:
    sys.exit(
        "ERROR: vtk package not found.\n"
        "Install with:  pip install vtk"
    )


_LJTS_SAT: dict[float, tuple[float, float]] = {
    0.3: (0.840, 0.0004),
    0.4: (0.815, 0.0020),
    0.6: (0.752, 0.0120),
    0.7: (0.714, 0.0230),
    0.8: (0.668, 0.0430),
    0.9: (0.608, 0.0820),
    1.0: (0.530, 0.1500),
}


def get_interface_density(T: float) -> float:
    """
    Return  rho_interface = (rho_liq + rho_vap) / 2  for temperature T.
    Uses linear interpolation between tabulated LJTS values.
    """
    temps = sorted(_LJTS_SAT.keys())
    T = float(np.clip(T, temps[0], temps[-1]))
    for i in range(len(temps) - 1):
        if temps[i] <= T <= temps[i + 1]:
            f = (T - temps[i]) / (temps[i + 1] - temps[i])
            rl = _LJTS_SAT[temps[i]][0] * (1 - f) + _LJTS_SAT[temps[i + 1]][0] * f
            rv = _LJTS_SAT[temps[i]][1] * (1 - f) + _LJTS_SAT[temps[i + 1]][1] * f
            return (rl + rv) / 2.0
    rl, rv = _LJTS_SAT[temps[0]]
    return (rl + rv) / 2.0




def load_data(path: str) -> tuple[np.ndarray, Optional[np.ndarray]]:
  
    path = str(path)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"File not found: {path}")

    if path.lower().endswith(".pvtu"):
        reader = vtk.vtkXMLPUnstructuredGridReader()
    else:
        reader = vtk.vtkXMLUnstructuredGridReader()

    reader.SetFileName(path)
    reader.Update()
    output = reader.GetOutput()

    if output is None:
        raise ValueError(f"VTK reader returned None for {path}")

    pts = output.GetPoints()
    if pts is None or pts.GetNumberOfPoints() == 0:
        raise ValueError(f"No point data in {path}")

    positions = vtk_to_numpy(pts.GetData()).astype(np.float64).copy()

    type_array = output.GetPointData().GetArray("typeIds")
    if type_array is not None:
        type_ids = vtk_to_numpy(type_array).astype(np.int32).copy()
    else:
        warnings.warn(
            f"'typeIds' point array absent in {os.path.basename(path)}. "
            "Will fall back to z-coordinate filtering.",
            stacklevel=2,
        )
        type_ids = None

    return positions, type_ids


def filter_droplet_particles(
    positions: np.ndarray,
    type_ids: Optional[np.ndarray],
    droplet_type: int,
    surface_z: float,
    fallback_margin: float = 1.0,
) -> np.ndarray:
   
    if len(positions) == 0:
        raise ValueError("positions array is empty.")

    if type_ids is not None:
        mask = type_ids == droplet_type
        n_matched = int(mask.sum())
        if n_matched > 0:
            return positions[mask]
        warnings.warn(
            f"No particles found with typeId={droplet_type}. "
            "Falling back to z-coordinate filter.",
            stacklevel=2,
        )

    # z-coordinate fallback
    z_threshold = surface_z + fallback_margin
    mask = positions[:, 2] > z_threshold
    if mask.sum() == 0:
        raise ValueError(
            f"Fallback z-filter (z > {z_threshold:.2f}) returned no particles. "
            "Check --surface-z."
        )
    return positions[mask]

def find_droplet_center(positions: np.ndarray) -> tuple[float, float]:
   
    return float(np.mean(positions[:, 0])), float(np.mean(positions[:, 1]))


def extract_central_slice(
    positions: np.ndarray,
    center_y: float,
    half_width: float,
) -> np.ndarray:
   
    mask = np.abs(positions[:, 1] - center_y) < half_width
    sliced = positions[mask]
    if len(sliced) < 10:
        raise ValueError(
            f"Only {len(sliced)} particles in the central slab "
            f"(|y - {center_y:.1f}| < {half_width:.1f} σ). "
            "Increase --slice-width."
        )
    return sliced

def build_2d_density(
    slice_positions: np.ndarray,
    surface_z: float,
    bin_size: float = 0.5,
    smooth_sigma: float = 1.0,
    slice_half_width: float = 5.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
   
    x = slice_positions[:, 0]
    z = slice_positions[:, 2] - surface_z   # shift so surface = 0

    above = z > 0.0
    x, z = x[above], z[above]
    if len(x) == 0:
        raise ValueError("No slice particles found above the substrate (z > 0).")

    x_min, x_max = x.min() - bin_size, x.max() + bin_size
    z_max = z.max() + bin_size

    x_bins = np.arange(x_min, x_max + bin_size, bin_size)
    z_bins = np.arange(0.0, z_max + bin_size, bin_size)

    counts, _, _ = np.histogram2d(x, z, bins=[x_bins, z_bins])

    cell_volume = bin_size * bin_size * 2.0 * slice_half_width
    density = counts / cell_volume

    density = gaussian_filter(density, sigma=smooth_sigma)
    return density, x_bins, z_bins



def _remove_outliers_iqr(
    x: np.ndarray,
    z: np.ndarray,
    factor: float = 1.5,
) -> tuple[np.ndarray, np.ndarray]:
   
    if len(x) < 4:
        return x, z

    z_lo_thresh = z.min() + 0.10 * (z.max() - z.min())
    upper_mask = z > z_lo_thresh

    if upper_mask.sum() < 4:
        return x, z

    q1, q3 = np.percentile(x[upper_mask], [25, 75])
    iqr = q3 - q1
    if iqr < 1e-10:
        return x, z

    # Keep points in the lower zone unconditionally; filter the upper zone.
    keep = (~upper_mask) | ((x >= q1 - factor * iqr) & (x <= q3 + factor * iqr))
    return x[keep], z[keep]


def detect_interface(
    density: np.ndarray,
    x_bins: np.ndarray,
    z_bins: np.ndarray,
    rho_interface: float,
    z_fit_min: float = 2.0,
    outlier_iqr_factor: float = 1.5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    
    x_mid = 0.5 * (x_bins[:-1] + x_bins[1:])
    z_mid = 0.5 * (z_bins[:-1] + z_bins[1:])
    dx = float(x_bins[1] - x_bins[0])

    left_x_list:  list[float] = []
    left_z_list:  list[float] = []
    right_x_list: list[float] = []
    right_z_list: list[float] = []

    for iz, z_val in enumerate(z_mid):
        if z_val < z_fit_min:
            continue

        row = density[:, iz]
        above = np.where(row >= rho_interface)[0]
        if len(above) == 0:
            continue

        left_idx  = int(above[0])
        right_idx = int(above[-1])

        # ── left boundary: interpolate between (left_idx-1) and left_idx ──
        if left_idx > 0:
            rho_prev = row[left_idx - 1]
            rho_curr = row[left_idx]
            denom = rho_curr - rho_prev
            if abs(denom) > 1e-30:
                f = (rho_interface - rho_prev) / denom
                lx = x_mid[left_idx - 1] + f * dx
            else:
                lx = x_mid[left_idx]
        else:
            lx = x_mid[left_idx]

        # ── right boundary: interpolate between right_idx and (right_idx+1) ──
        if right_idx + 1 < len(x_mid):
            rho_curr = row[right_idx]
            rho_next = row[right_idx + 1]
            denom = rho_curr - rho_next          # positive: density drops going right
            if abs(denom) > 1e-30:
                f = (rho_curr - rho_interface) / denom
                rx = x_mid[right_idx] + f * dx
            else:
                rx = x_mid[right_idx]
        else:
            rx = x_mid[right_idx]

        left_x_list.append(lx)
        left_z_list.append(z_val)
        right_x_list.append(rx)
        right_z_list.append(z_val)

    left_x  = np.array(left_x_list,  dtype=np.float64)
    left_z  = np.array(left_z_list,  dtype=np.float64)
    right_x = np.array(right_x_list, dtype=np.float64)
    right_z = np.array(right_z_list, dtype=np.float64)

    left_x,  left_z  = _remove_outliers_iqr(left_x,  left_z,  outlier_iqr_factor)
    right_x, right_z = _remove_outliers_iqr(right_x, right_z, outlier_iqr_factor)

    return left_x, left_z, right_x, right_z


# ─────────────────────────────────────────────────────────────────────────────
# 7.  CIRCLE FIT
# ─────────────────────────────────────────────────────────────────────────────

def fit_circle(
    x: np.ndarray,
    z: np.ndarray,
) -> tuple[float, float, float, float]:
    
    if len(x) < 3:
        raise ValueError(f"Need ≥ 3 points to fit a circle; got {len(x)}.")

    # ── algebraic initial estimate ──────────────────────────────────────────
    # Solve  [x  z  1] · [a  b  c]ᵀ = −(x² + z²)  in the least-squares sense.
    A = np.column_stack([x, z, np.ones(len(x))])
    b_vec = -(x ** 2 + z ** 2)
    try:
        coeffs, _, _, _ = np.linalg.lstsq(A, b_vec, rcond=None)
        a_c, b_c, c_c = coeffs
        cx0 = -a_c / 2.0
        cz0 = -b_c / 2.0
        r0_sq = cx0 ** 2 + cz0 ** 2 - c_c
        R0 = float(np.sqrt(max(r0_sq, 1e-6)))
    except np.linalg.LinAlgError:
        cx0 = float(np.mean(x))
        cz0 = float(np.mean(z))
        R0 = float(np.std(x) + np.std(z))

    # ── geometric Levenberg-Marquardt refinement ────────────────────────────
    def _residuals(p: np.ndarray) -> np.ndarray:
        cx_, cz_, R_ = p
        return np.sqrt((x - cx_) ** 2 + (z - cz_) ** 2) - R_

    try:
        result = least_squares(
            _residuals,
            x0=[cx0, cz0, R0],
            method="lm",
            max_nfev=20_000,
            ftol=1e-10,
            xtol=1e-10,
        )
        cx, cz, R = float(result.x[0]), float(result.x[1]), float(result.x[2])
        rms = float(np.sqrt(np.mean(result.fun ** 2)))
    except Exception as exc:
        warnings.warn(f"Geometric circle-fit refinement failed ({exc}); using algebraic solution.", stacklevel=2)
        cx, cz, R = cx0, cz0, R0
        rms = float(np.sqrt(np.mean(_residuals(np.array([cx, cz, R])) ** 2)))

    R = abs(R)
    if R < 1e-6:
        raise ValueError(f"Circle fit returned near-zero radius R={R:.6f}.")

    return cx, cz, R, rms


def compute_contact_angles(
    cx: float,
    cz: float,
    R: float,
) -> tuple[float, float, float, float, float]:
    
    if R <= 0.0:
        raise ValueError(f"Invalid radius R = {R:.4f}")
    if abs(cz) > R:
        # The fitted circle does not geometrically intersect z = 0.
        # This occurs for highly hydrophobic droplets (θ → 180°) where the
        # droplet barely touches the surface.  Clamp cos(θ) to [-1, 1] and
        # place the contact points at cx (zero lateral spread).
        warnings.warn(
            f"Circle does not intersect z = 0 (|cz|={abs(cz):.3f} > R={R:.3f}). "
            f"Droplet is super-hydrophobic — clamping θ to {np.degrees(np.arccos(np.clip(-cz/R,-1,1))):.1f}°.",
            stacklevel=3,
        )
        discriminant = 0.0
    else:
        discriminant = float(np.sqrt(max(R ** 2 - cz ** 2, 0.0)))
    x_left_cp  = cx - discriminant
    x_right_cp = cx + discriminant

    # Both left and right angles are mathematically identical for a single
    # circle fit (the formula depends only on cz and R, not on the cx offset).
    # They differ only when the left and right interface halves are fitted
    # separately, which is done in process_file().
    cos_theta = float(np.clip(-cz / R, -1.0, 1.0))
    theta = float(np.degrees(np.arccos(cos_theta)))

    return theta, theta, theta, float(x_left_cp), float(x_right_cp)


def _draw_contact_angle_arc(
    ax: plt.Axes,
    x_cp: float,
    z_cp: float,
    theta: float,
    R: float,
    color: str,
    side: str,
    label: str,
) -> None:
    
    arc_r = max(R * 0.10, 2.5)

    if side == "left":
        theta1, theta2 = 0.0, theta
        mid_deg = theta / 2.0
        label_ha = "right"
    else:
        theta1 = 180.0 - theta
        theta2 = 180.0
        mid_deg = 180.0 - theta / 2.0
        label_ha = "left"

    arc = mpatches.Arc(
        (x_cp, z_cp),
        2.0 * arc_r,
        2.0 * arc_r,
        angle=0.0,
        theta1=theta1,
        theta2=theta2,
        color=color,
        lw=2.5,
        zorder=12,
    )
    ax.add_patch(arc)

    mid_rad = np.radians(mid_deg)
    lx = x_cp + arc_r * 2.1 * np.cos(mid_rad)
    lz = z_cp + arc_r * 2.1 * np.sin(mid_rad)
    ax.annotate(
        label,
        xy=(lx, lz),
        fontsize=12,
        color=color,
        fontweight="bold",
        ha=label_ha,
        va="center",
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=color, alpha=0.90),
        zorder=13,
    )


def create_diagnostic_plot(
    slice_positions: np.ndarray,
    density: np.ndarray,
    x_bins: np.ndarray,
    z_bins: np.ndarray,
    left_x: np.ndarray,
    left_z: np.ndarray,
    right_x: np.ndarray,
    right_z: np.ndarray,
    cx: float,
    cz: float,
    R: float,
    theta_left: float,
    theta_right: float,
    theta_avg: float,
    x_left_cp: float,
    x_right_cp: float,
    surface_z: float,
    z_fit_min: float,
    title: str,
    out_path: str,
    dpi: int = 300,
) -> None:
   
    fig, ax = plt.subplots(figsize=(11, 9))
    fig.patch.set_facecolor("white")

    # ── density heatmap ──────────────────────────────────────────────────────
    Xg, Zg = np.meshgrid(x_bins, z_bins)
    D = density.T          # shape (nz, nx) for pcolormesh
    pos_vals = density[density > 0]
    vmax = float(np.percentile(pos_vals, 95)) if len(pos_vals) > 0 else 1.0
    im = ax.pcolormesh(
        Xg, Zg, D,
        cmap="Blues", vmin=0.0, vmax=vmax,
        shading="flat", alpha=0.70, rasterized=True, zorder=1,
    )
    cbar = plt.colorbar(im, ax=ax, label=r"$\rho\ /\ \sigma^{-3}$",
                        shrink=0.60, pad=0.02)
    cbar.ax.tick_params(labelsize=11)

    # ── substrate ────────────────────────────────────────────────────────────
    x_lo, x_hi = float(x_bins[0]), float(x_bins[-1])
    z_top = float(z_bins[-1])
    ax.fill_between([x_lo, x_hi], -3.5, 0.0,
                    color="#8B6914", alpha=0.55, zorder=2, label="Substrate")
    ax.axhline(0.0, color="#5C4000", lw=2.0, zorder=3)

    # ── wall-exclusion zone ──────────────────────────────────────────────────
    if z_fit_min > 0:
        ax.axhspan(0.0, z_fit_min, color="lightsalmon", alpha=0.28, zorder=1,
                   label=f"Excl. zone  z < {z_fit_min:.1f} σ")

    # ── particle scatter (central slice) ─────────────────────────────────────
    z_shifted = slice_positions[:, 2] - surface_z
    above = z_shifted > 0.0
    ax.scatter(
        slice_positions[above, 0], z_shifted[above],
        s=2.0, c="steelblue", alpha=0.20, linewidths=0,
        rasterized=True, zorder=4, label="Slice particles",
    )

    # ── interface points ─────────────────────────────────────────────────────
    ax.scatter(left_x,  left_z,  s=22, c="royalblue",   zorder=8,
               edgecolors="navy",        linewidths=0.5, label="Interface (left)")
    ax.scatter(right_x, right_z, s=22, c="darkorange",  zorder=8,
               edgecolors="saddlebrown", linewidths=0.5, label="Interface (right)")

    # ── fitted circle arc (z ≥ 0 only) ───────────────────────────────────────
    phi = np.linspace(0.0, 2.0 * np.pi, 1800)
    xc_all = cx + R * np.cos(phi)
    zc_all = cz + R * np.sin(phi)
    valid  = zc_all >= 0.0
    xc_plot = np.where(valid, xc_all, np.nan)
    zc_plot = np.where(valid, zc_all, np.nan)
    ax.plot(xc_plot, zc_plot, color="crimson", lw=2.8, zorder=9,
            label=f"Circle fit  R = {R:.2f} σ")
    ax.plot(cx, cz, "P", color="crimson", ms=10, zorder=10,
            label=f"Centre  ({cx:.1f}, {cz:.1f}) σ")

    # ── contact points ────────────────────────────────────────────────────────
    ax.plot(x_left_cp,  0.0, "v", color="darkgreen",  ms=11, zorder=11,
            label=f"Left  contact  x = {x_left_cp:.1f} σ",
            markeredgecolor="white", markeredgewidth=0.5)
    ax.plot(x_right_cp, 0.0, "v", color="darkorchid", ms=11, zorder=11,
            label=f"Right contact  x = {x_right_cp:.1f} σ",
            markeredgecolor="white", markeredgewidth=0.5)

    # ── contact-angle arcs ────────────────────────────────────────────────────
    _draw_contact_angle_arc(ax, x_left_cp,  0.0, theta_left,  R,
                            "darkgreen",  "left",
                            f"θ_L = {theta_left:.1f}°")
    _draw_contact_angle_arc(ax, x_right_cp, 0.0, theta_right, R,
                            "darkorchid", "right",
                            f"θ_R = {theta_right:.1f}°")

    # ── summary annotation ────────────────────────────────────────────────────
    wetting = "hydrophilic" if theta_avg < 90.0 else (
              "superhydrophobic" if theta_avg >= 150.0 else "hydrophobic")
    ax.text(
        cx, z_top * 0.87,
        (f"θ_avg = {theta_avg:.1f}°  ({wetting})\n"
         f"θ_L = {theta_left:.1f}°     θ_R = {theta_right:.1f}°"),
        ha="center", va="center", fontsize=13, fontweight="bold",
        color="white",
        bbox=dict(boxstyle="round,pad=0.45", fc="navy", ec="none", alpha=0.82),
        zorder=14,
    )

    # ── formatting ────────────────────────────────────────────────────────────
    ax.set_xlabel(r"$x\ /\ \sigma$", fontsize=15)
    ax.set_ylabel(r"$z\ /\ \sigma$", fontsize=15)
    ax.set_title(title, fontsize=13, pad=10)
    ax.set_xlim(x_bins[0], x_bins[-1])
    ax.set_ylim(-3.5, z_top)
    ax.set_aspect("equal")
    ax.tick_params(labelsize=12)
    ax.grid(True, ls="--", lw=0.5, alpha=0.40, color="gray", zorder=0)
    ax.legend(
        fontsize=9, loc="upper right", framealpha=0.92,
        edgecolor="gray", ncol=2,
    )

    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def save_results(results: list[dict], out_dir: str) -> None:
    
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, "contact_angles.csv")
    txt_path = os.path.join(out_dir, "contact_angles.txt")

    fieldnames = [
        "filename", "timestep",
        "left_angle", "right_angle", "average_angle",
        "radius", "center_x", "center_z", "fit_rms",
    ]

    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow({k: r.get(k, "") for k in fieldnames})

    col_w = 28
    sep = "=" * (col_w + 10 * 9)
    with open(txt_path, "w") as fh:
        fh.write(f"Static Contact Angle Results\n{sep}\n")
        fh.write(
            f"{'Filename':<{col_w}}  {'Step':>6}  "
            f"{'θ_L':>7}  {'θ_R':>7}  {'θ_avg':>7}  "
            f"{'R':>7}  {'cx':>8}  {'cz':>8}  {'RMS':>7}\n"
        )
        fh.write("-" * (col_w + 10 * 9) + "\n")
        for r in results:
            fh.write(
                f"{r['filename']:<{col_w}}  {r['timestep']:>6}  "
                f"{r['left_angle']:>7.2f}  {r['right_angle']:>7.2f}  "
                f"{r['average_angle']:>7.2f}  "
                f"{r['radius']:>7.2f}  {r['center_x']:>8.3f}  "
                f"{r['center_z']:>8.3f}  {r['fit_rms']:>7.4f}\n"
            )
        fh.write(f"{sep}\n")
        if results:
            avgs = [r["average_angle"] for r in results]
            fh.write(
                f"\nMean θ_avg = {np.mean(avgs):.2f}°  "
                f"(std = {np.std(avgs):.2f}°,  "
                f"min = {np.min(avgs):.2f}°,  "
                f"max = {np.max(avgs):.2f}°)\n"
            )

    print(f"  CSV     → {csv_path}")
    print(f"  Summary → {txt_path}")


def process_file(vtu_path: str, args: argparse.Namespace, out_dir: str,
                 file_index: int) -> Optional[dict]:
    """
    Run the full contact-angle measurement pipeline for a single snapshot.

    Returns a dict of results, or None if the file cannot be processed.
    """
    fname = os.path.basename(vtu_path)
    print(f"\n{'─'*68}")
    print(f"  [{file_index:04d}]  {fname}")
    print(f"{'─'*68}")

    # ── infer temperature from folder/filename (e.g. temp07 → T=0.7) ────────
    T = args.temperature
    m = re.search(r"temp(\d+)", vtu_path)
    if m:
        T = int(m.group(1)) / 10.0

    rho_interface = get_interface_density(T)
    print(f"  T = {T:.2f}   rho_interface = {rho_interface:.5f} σ⁻³")

    # ── extract integer timestep from filename ───────────────────────────────
    ts_match = re.search(r"_(\d+)\.p?vtu$", fname)
    timestep = int(ts_match.group(1)) if ts_match else file_index

    # ── step 1: load data ────────────────────────────────────────────────────
    try:
        positions, type_ids = load_data(vtu_path)
    except (FileNotFoundError, ValueError) as exc:
        print(f"  [SKIP] Load failed: {exc}")
        return None
    print(f"  Total particles   : {len(positions):>8,}")

    # ── step 2: filter droplet particles ────────────────────────────────────
    try:
        droplet = filter_droplet_particles(
            positions, type_ids,
            args.droplet_type, args.surface_z,
        )
    except ValueError as exc:
        print(f"  [SKIP] Particle filter failed: {exc}")
        return None
    print(f"  Droplet particles : {len(droplet):>8,}")

    # ── step 3: droplet centre ───────────────────────────────────────────────
    cx0, cy0 = find_droplet_center(droplet)
    print(f"  Droplet centroid  : x = {cx0:.2f} σ,  y = {cy0:.2f} σ")

    # ── step 4: central slice ────────────────────────────────────────────────
    try:
        sliced = extract_central_slice(droplet, cy0, args.slice_half_width)
    except ValueError as exc:
        print(f"  [SKIP] Slice failed: {exc}")
        return None
    print(f"  Slice particles   : {len(sliced):>8,}  "
          f"(±{args.slice_half_width:.1f} σ around y = {cy0:.1f} σ)")

    # ── step 5: 2-D density histogram ────────────────────────────────────────
    try:
        density, x_bins, z_bins = build_2d_density(
            sliced, args.surface_z,
            bin_size=args.bin_size,
            smooth_sigma=args.smooth_sigma,
            slice_half_width=args.slice_half_width,
        )
    except ValueError as exc:
        print(f"  [SKIP] Density build failed: {exc}")
        return None

    # ── step 6: interface detection ──────────────────────────────────────────
    left_x, left_z, right_x, right_z = detect_interface(
        density, x_bins, z_bins,
        rho_interface=rho_interface,
        z_fit_min=args.z_fit_min,
    )
    n_iface = len(left_x) + len(right_x)
    print(f"  Interface points  : {n_iface}  "
          f"(left: {len(left_x)}, right: {len(right_x)})")

    if len(left_x) < args.min_interface_pts or len(right_x) < args.min_interface_pts:
        print(f"  [SKIP] Too few interface points (min required: {args.min_interface_pts}).")
        return None

    # ── step 7: circle fits ──────────────────────────────────────────────────
    # Global fit: all interface points combined.
    all_x = np.concatenate([left_x, right_x])
    all_z = np.concatenate([left_z, right_z])
    try:
        cx_g, cz_g, R_g, rms_g = fit_circle(all_x, all_z)
        print(f"  Global fit: cx={cx_g:.3f}  cz={cz_g:.3f}  R={R_g:.3f}  RMS={rms_g:.4f}")
    except ValueError as exc:
        print(f"  [SKIP] Global circle fit failed: {exc}")
        return None

    # Split fits: left half → theta_L,  right half → theta_R.
    # This captures genuine asymmetry on structured surfaces.
    try:
        _, cz_L, R_L, rms_L = fit_circle(left_x, left_z)
        cos_L = float(np.clip(-cz_L / R_L, -1.0, 1.0))
        theta_L = float(np.degrees(np.arccos(cos_L)))
        print(f"  Left  fit : cz={cz_L:.3f}  R={R_L:.3f}  θ_L={theta_L:.2f}°")
    except (ValueError, ZeroDivisionError) as exc:
        warnings.warn(f"Left half fit failed ({exc}); using global fit angle.", stacklevel=2)
        theta_L = None
        rms_L   = rms_g

    try:
        _, cz_R, R_R, rms_R = fit_circle(right_x, right_z)
        cos_R = float(np.clip(-cz_R / R_R, -1.0, 1.0))
        theta_R = float(np.degrees(np.arccos(cos_R)))
        print(f"  Right fit : cz={cz_R:.3f}  R={R_R:.3f}  θ_R={theta_R:.2f}°")
    except (ValueError, ZeroDivisionError) as exc:
        warnings.warn(f"Right half fit failed ({exc}); using global fit angle.", stacklevel=2)
        theta_R = None
        rms_R   = rms_g

    # ── step 8: contact angles ───────────────────────────────────────────────
    try:
        theta_g, _, _, x_left_cp, x_right_cp = compute_contact_angles(cx_g, cz_g, R_g)
    except ValueError as exc:
        print(f"  [SKIP] Contact angle failed: {exc}")
        return None

    # Use split-fit angles when available, global otherwise.
    theta_left  = theta_L if theta_L is not None else theta_g
    theta_right = theta_R if theta_R is not None else theta_g
    theta_avg   = (theta_left + theta_right) / 2.0

    wetting = ("hydrophilic"      if theta_avg < 90.0 else
               "superhydrophobic" if theta_avg >= 150.0 else
               "hydrophobic")
    print(f"\n  ★  θ_L = {theta_left:.2f}°   θ_R = {theta_right:.2f}°   "
          f"θ_avg = {theta_avg:.2f}°   [{wetting}]")

    # ── step 9: diagnostic plot ──────────────────────────────────────────────
    png_name = f"contact_angle_{file_index:04d}.png"
    out_path = os.path.join(out_dir, png_name)
    title = (f"{fname}  |  timestep {timestep}  |  "
             f"θ_avg = {theta_avg:.1f}°  ({wetting})")
    try:
        create_diagnostic_plot(
            slice_positions=sliced,
            density=density,
            x_bins=x_bins,
            z_bins=z_bins,
            left_x=left_x,
            left_z=left_z,
            right_x=right_x,
            right_z=right_z,
            cx=cx_g,
            cz=cz_g,
            R=R_g,
            theta_left=theta_left,
            theta_right=theta_right,
            theta_avg=theta_avg,
            x_left_cp=x_left_cp,
            x_right_cp=x_right_cp,
            surface_z=args.surface_z,
            z_fit_min=args.z_fit_min,
            title=title,
            out_path=out_path,
            dpi=args.dpi,
        )
        print(f"  Figure → {out_path}")
    except Exception as exc:
        warnings.warn(f"Figure generation failed: {exc}", stacklevel=2)

    return {
        "filename":      fname,
        "timestep":      timestep,
        "left_angle":    round(theta_left,  4),
        "right_angle":   round(theta_right, 4),
        "average_angle": round(theta_avg,   4),
        "radius":        round(R_g,         4),
        "center_x":      round(cx_g,        4),
        "center_z":      round(cz_g,        4),
        "fit_rms":       round(rms_g,       6),
    }

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Static contact-angle analysis for sessile MD droplet snapshots.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # ── input ─────────────────────────────────────────────────────────────────
    p.add_argument(
        "--input", "-i",
        default="/Users/melisaaslan/Desktop",
        metavar="DIR",
        help="Root directory to search for VTU/PVTU files.",
    )
    p.add_argument(
        "--pattern", "-p",
        default="**/*_Particles_*.pvtu",
        metavar="GLOB",
        help="Glob pattern (relative to --input) to select snapshot files. "
             "Default targets .pvtu master files (parallel runs). "
             "For single-process runs use '**/*_Particles_*.vtu'.",
    )
    p.add_argument(
        "--last-n",
        type=int,
        default=10,
        metavar="N",
        help="Average the contact angle over the last N timesteps per simulation run. "
             "Use 1 to analyse only the final snapshot.",
    )
    p.add_argument(
        "--output", "-o",
        default=None,
        metavar="DIR",
        help="Output directory for PNGs, CSV, and summary text. "
             "Defaults to <input>/contact_angle_results/.",
    )

    # ── physics ───────────────────────────────────────────────────────────────
    p.add_argument(
        "--surface-z",
        type=float,
        default=13.0,
        metavar="SIGMA",
        help="z-coordinate of the top of the solid substrate (sigma). "
             "Increase for boss/pit surfaces.",
    )
    p.add_argument(
        "--temperature", "-T",
        type=float,
        default=0.7,
        metavar="T*",
        help="Reduced LJTS temperature for the Gibbs dividing-surface density. "
             "Overridden automatically when 'tempNN' appears in the folder name.",
    )
    p.add_argument(
        "--droplet-type",
        type=int,
        default=0,
        metavar="INT",
        help="Integer typeId that identifies liquid-droplet particles.",
    )

    # ── binning ───────────────────────────────────────────────────────────────
    p.add_argument(
        "--bin-size",
        type=float,
        default=0.5,
        metavar="SIGMA",
        help="Histogram cell size (sigma) used for the 2-D density field.",
    )
    p.add_argument(
        "--smooth-sigma",
        type=float,
        default=1.0,
        metavar="BINS",
        help="Gaussian smoothing kernel width (in bin units) for the density.",
    )
    p.add_argument(
        "--slice-half-width",
        type=float,
        default=10.0,
        metavar="SIGMA",
        help="Half-thickness (sigma) of the central y-slab used for the 2-D projection.",
    )
    p.add_argument(
        "--z-fit-min",
        type=float,
        default=2.0,
        metavar="SIGMA",
        help="Minimum z (sigma) above the substrate used for circle fitting. "
             "Excludes wall-layering effects (Becker et al. 2014).",
    )
    p.add_argument(
        "--min-interface-pts",
        type=int,
        default=8,
        metavar="N",
        help="Minimum number of detected interface points per side required "
             "for a reliable circle fit.",
    )

    # ── output ────────────────────────────────────────────────────────────────
    p.add_argument(
        "--dpi",
        type=int,
        default=300,
        metavar="INT",
        help="Output figure resolution in DPI.",
    )

    return p


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    # ── resolve paths ─────────────────────────────────────────────────────────
    input_dir  = os.path.abspath(args.input)
    output_dir = args.output or os.path.join(input_dir, "contact_angle_results")
    os.makedirs(output_dir, exist_ok=True)

    # ── discover files ────────────────────────────────────────────────────────
    search_pattern = os.path.join(input_dir, args.pattern)
    all_files = sorted(glob.glob(search_pattern, recursive=True))

    # Exclude piece files inside data/ subdirectories — these are per-rank fragments
    # produced by parallel runs and only contain a fraction of particles each.
    # The master .pvtu files at the top level reference all pieces and are correct.
    all_files = [f for f in all_files
                 if os.path.basename(os.path.dirname(f)) != "data"]

    if not all_files:
        print(f"No files found matching: {search_pattern}")
        print("  Tip: for parallel runs use --pattern '**/*_Particles_*.pvtu'")
        sys.exit(1)

    # Keep the last N timesteps per simulation run.
    # Filename pattern: <prefix>_Particles_<NNNNNN>.<ext>
    # Group by (directory + prefix), sort by timestep, take last --last-n.
    _ts_re = re.compile(r"^(.*_Particles_)(\d+)(\.[^.]+)$", re.IGNORECASE)
    _groups: dict[str, list[tuple[int, str]]] = {}  # key -> [(timestep, path), ...]
    _ungrouped: list[str] = []
    for path in all_files:
        basename = os.path.basename(path)
        m = _ts_re.match(basename)
        if m:
            key = os.path.join(os.path.dirname(path), m.group(1))
            ts = int(m.group(2))
            _groups.setdefault(key, []).append((ts, path))
        else:
            _ungrouped.append(path)

    # For each simulation, sort by timestep and keep only the last N.
    last_n = args.last_n
    _sim_files: dict[str, list[str]] = {}  # key -> [path, ...] (sorted, last N)
    for key, entries in _groups.items():
        entries.sort(key=lambda x: x[0])
        _sim_files[key] = [p for _, p in entries[-last_n:]]

    files = [p for paths in _sim_files.values() for p in paths] + sorted(_ungrouped)

    n_sims = len(_sim_files)
    print(f"\nFound {len(all_files)} file(s) total; "
          f"keeping last {last_n} timestep(s) × {n_sims} simulation(s) = {len(files)} file(s).")
    print(f"Output directory: {output_dir}\n")

    # ── process each file ─────────────────────────────────────────────────────
    all_results: list[dict] = []
    skipped: list[tuple[str, str]] = []  # (filename, reason)
    for idx, vtu_path in enumerate(files, start=1):
        result = process_file(vtu_path, args, output_dir, idx)
        if result is not None:
            all_results.append(result)
        else:
            skipped.append((os.path.basename(vtu_path), "see output above"))

    # ── save CSV + summary ────────────────────────────────────────────────────
    if all_results:
        save_results(all_results, output_dir)
    else:
        print("\n[WARNING] No files were successfully processed.")

    # ── console summary ───────────────────────────────────────────────────────
    W = 50
    print(f"\n{'═'*68}")
    print("  FINAL SUMMARY  –  Static Contact Angle Results")
    print(f"{'═'*68}")
    print(f"  Files found   : {len(files)}")
    print(f"  Processed OK  : {len(all_results)}")
    print(f"  Skipped       : {len(skipped)}")
    print()
    if all_results:
        # Per-timestep table
        print(f"  {'Filename':<{W}}  {'θ_L':>7}  {'θ_R':>7}  {'θ_avg':>7}")
        print(f"  {'─'*W}  {'─'*7}  {'─'*7}  {'─'*7}")
        for r in all_results:
            fname = r['filename']
            if len(fname) > W:
                fname = "…" + fname[-(W-1):]
            print(f"  {fname:<{W}}  "
                  f"{r['left_angle']:>7.2f}  "
                  f"{r['right_angle']:>7.2f}  "
                  f"{r['average_angle']:>7.2f}")

        # Per-simulation average over last N timesteps
        print(f"\n  {'─'*68}")
        print(f"  AVERAGED OVER LAST {last_n} TIMESTEP(S) PER SIMULATION")
        print(f"  {'─'*68}")
        # Group results back by simulation key
        _res_by_sim: dict[str, list[dict]] = {}
        for r in all_results:
            # Match result filename back to a simulation key
            for key, paths in _sim_files.items():
                if r['filename'] in [os.path.basename(p) for p in paths]:
                    _res_by_sim.setdefault(key, []).append(r)
                    break
        for key, res_list in _res_by_sim.items():
            sim_name = os.path.basename(key).rstrip("_")
            if len(sim_name) > W:
                sim_name = "…" + sim_name[-(W-1):]
            tl = np.mean([r["left_angle"]    for r in res_list])
            tr = np.mean([r["right_angle"]   for r in res_list])
            ta = np.mean([r["average_angle"] for r in res_list])
            st = np.std( [r["average_angle"] for r in res_list])
            print(f"  {sim_name:<{W}}  "
                  f"{tl:>7.2f}  {tr:>7.2f}  {ta:>7.2f}  (std={st:.2f}°, n={len(res_list)})")
    if skipped:
        print(f"\n  {'─'*68}")
        print("  SKIPPED FILES  (check the per-file output above for details)")
        print(f"  {'─'*68}")
        for fname, reason in skipped:
            short = fname if len(fname) <= W else "…" + fname[-(W-1):]
            print(f"  [SKIP]  {short}")
    print(f"{'═'*68}\n")


if __name__ == "__main__":
    main()
