"""
contact_angle_paraview.py
=========================
Static contact angle measurement — runs inside ParaView's Python Shell.

USAGE
-----
1. Load your PVTU/VTU file in ParaView and make it the active source.
2. Open  View → Python Shell
3. Click  "Run Script"  and select this file.

SETTINGS — edit the block below before running.
"""

# ── USER SETTINGS ─────────────────────────────────────────────────────────────
SURFACE_Z       = 12.349   # z-coordinate of the top of the surface (sigma)
                           # Boss:  12.349  |  Substrate only:  4.266
TEMPERATURE     = None     # Set e.g. 0.6, or None to auto-detect from filename
DROPLET_TYPE_ID = 0        # typeId of droplet particles (0 = fluid)
SLICE_HALF_WIDTH = 5.0     # half-width of the central y-slab (sigma)
BIN_SIZE        = 0.5      # density histogram bin size (sigma)
SMOOTH_SIGMA    = 1.5      # Gaussian smoothing (bins)
Z_FIT_MIN       = 2.0      # exclude interface points below this z (sigma)
MIN_IFACE_PTS   = 8        # minimum interface points per side
# ──────────────────────────────────────────────────────────────────────────────

import re
import warnings
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.optimize import least_squares

try:
    from paraview.simple import GetActiveSource
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    _IN_PARAVIEW = True
except ImportError:
    _IN_PARAVIEW = False

# ── LJTS saturation densities (Vrabec et al. Mol. Phys. 2006) ─────────────────
_LJTS_SAT = {
    0.3: (0.840, 0.0004),
    0.4: (0.815, 0.0020),
    0.6: (0.752, 0.0120),
    0.7: (0.714, 0.0230),
    0.8: (0.668, 0.0430),
    0.9: (0.608, 0.0820),
    1.0: (0.530, 0.1500),
}


def _get_interface_density(T):
    temps = sorted(_LJTS_SAT.keys())
    T = float(np.clip(T, temps[0], temps[-1]))
    for i in range(len(temps) - 1):
        if temps[i] <= T <= temps[i + 1]:
            f = (T - temps[i]) / (temps[i + 1] - temps[i])
            rl = _LJTS_SAT[temps[i]][0] * (1-f) + _LJTS_SAT[temps[i+1]][0] * f
            rv = _LJTS_SAT[temps[i]][1] * (1-f) + _LJTS_SAT[temps[i+1]][1] * f
            return (rl + rv) / 2.0
    rl, rv = _LJTS_SAT[temps[0]]
    return (rl + rv) / 2.0


def _get_positions_and_types():
    """Read positions and typeIds from the active ParaView source."""
    src = GetActiveSource()
    if src is None:
        raise RuntimeError("No active source in ParaView. Load a file first.")

    src.UpdatePipeline()
    data = src.GetClientSideObject().GetOutput()

    # Handle composite datasets (PVTU → UnstructuredGrid pieces)
    if hasattr(data, 'GetBlock'):
        pieces = [data.GetBlock(i)
                  for i in range(data.GetNumberOfBlocks())
                  if data.GetBlock(i) is not None]
    else:
        pieces = [data]

    all_pts, all_types = [], []
    _type_array_name = None
    for piece in pieces:
        pts = vtk_to_numpy(piece.GetPoints().GetData())
        all_pts.append(pts)
        pd = piece.GetPointData()
        # Try common name variants
        arr = None
        for name in ("typeIds", "typeId", "type", "TypeId", "TypeIds"):
            arr = pd.GetArray(name)
            if arr is not None:
                _type_array_name = name
                break
        if arr is not None:
            all_types.append(vtk_to_numpy(arr).astype(int))
        else:
            all_types.append(np.full(len(pts), DROPLET_TYPE_ID, dtype=int))

    positions = np.vstack(all_pts)
    type_ids  = np.concatenate(all_types)
    return positions, type_ids


def _infer_temperature(src):
    """Try to read temperature from the source's filename."""
    if TEMPERATURE is not None:
        return float(TEMPERATURE)
    try:
        fname = src.GetProperty("FileName").GetElement(0)
    except Exception:
        try:
            fname = src.FileName
        except Exception:
            return 0.6  # fallback
    m = re.search(r"temp(\d+)", fname)
    if m:
        return int(m.group(1)) / 10.0
    return 0.6


def _fit_circle(x, z):
    """Algebraic init → Levenberg-Marquardt geometric refinement."""
    if len(x) < 3:
        raise ValueError(f"Need ≥3 points; got {len(x)}.")
    # Algebraic (Pratt) initial guess
    A = np.column_stack([2*x, 2*z, np.ones(len(x))])
    b = x**2 + z**2
    res, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    cx0, cz0 = res[0], res[1]
    R0 = float(np.sqrt(res[2] + cx0**2 + cz0**2))
    if R0 <= 0:
        R0 = 1.0

    def resid(p):
        return np.sqrt((x - p[0])**2 + (z - p[1])**2) - p[2]

    result = least_squares(resid, [cx0, cz0, R0], method="lm", max_nfev=2000)
    cx, cz, R = result.x
    if R <= 0:
        raise ValueError(f"Fit returned R={R:.4f}.")
    rms = float(np.sqrt(np.mean(resid(result.x)**2)))
    return float(cx), float(cz), float(R), rms


def _contact_angle(cz, R):
    """cos(θ) = -cz/R, clamped to [-1,1]."""
    cos_t = float(np.clip(-cz / R, -1.0, 1.0))
    if abs(cz) > R:
        warnings.warn(f"Super-hydrophobic: circle doesn't reach z=0 "
                      f"(|cz|={abs(cz):.2f} > R={R:.2f}). θ clamped.")
    return float(np.degrees(np.arccos(cos_t)))


def run():
    if not _IN_PARAVIEW:
        print("[ERROR] This script must be run inside ParaView's Python Shell.")
        return

    src = GetActiveSource()
    T = _infer_temperature(src)
    rho_iface = _get_interface_density(T)

    print(f"\n{'═'*60}")
    print("  Contact Angle Analysis  (ParaView)")
    print(f"{'═'*60}")
    print(f"  T = {T:.2f}   rho_interface = {rho_iface:.5f} σ⁻³")
    print(f"  surface_z = {SURFACE_Z:.3f} σ")

    # ── 1. load data ──────────────────────────────────────────────────────────
    positions, type_ids = _get_positions_and_types()
    print(f"  Total particles : {len(positions):,}")
    unique_ids, counts = np.unique(type_ids, return_counts=True)
    for uid, cnt in zip(unique_ids, counts):
        label = {0: "droplet", 1: "substrate", 2: "boss", 5: "boss"}.get(int(uid), f"typeId={uid}")
        print(f"    typeId={uid} ({label}): {cnt:,}")

    # ── 2. filter droplet particles ───────────────────────────────────────────
    mask = type_ids == DROPLET_TYPE_ID
    if mask.sum() == 0:
        print("  [WARN] No typeId==0 found; using z > surface_z + 1.0 as fallback.")
        mask = positions[:, 2] > (SURFACE_Z + 1.0)
    droplet = positions[mask]
    print(f"  Droplet particles: {len(droplet):,}")

    # ── 3. centroid ───────────────────────────────────────────────────────────
    cx0 = float(np.mean(droplet[:, 0]))
    cy0 = float(np.mean(droplet[:, 1]))
    print(f"  Centroid : x={cx0:.2f}  y={cy0:.2f} σ")

    # ── 4. central slab ───────────────────────────────────────────────────────
    slab_mask = np.abs(droplet[:, 1] - cy0) < SLICE_HALF_WIDTH
    slab = droplet[slab_mask]
    if len(slab) < 20:
        print(f"  [ERROR] Too few slab particles ({len(slab)}). "
              "Increase SLICE_HALF_WIDTH.")
        return
    print(f"  Slab particles : {len(slab):,}  (±{SLICE_HALF_WIDTH} σ)")

    # ── 5. 2D density histogram ───────────────────────────────────────────────
    x_arr = slab[:, 0]
    z_arr = slab[:, 2] - SURFACE_Z           # shift so surface = 0

    above = z_arr > 0.0
    x_arr, z_arr = x_arr[above], z_arr[above]

    x_bins = np.arange(x_arr.min(), x_arr.max() + BIN_SIZE, BIN_SIZE)
    z_bins = np.arange(0.0, z_arr.max() + BIN_SIZE, BIN_SIZE)
    density, _, _ = np.histogram2d(x_arr, z_arr,
                                   bins=[x_bins, z_bins])
    # normalise to number density (sigma^-3)
    vox = BIN_SIZE**2 * (2.0 * SLICE_HALF_WIDTH)
    density = density / vox
    density = gaussian_filter(density.T, sigma=SMOOTH_SIGMA)  # shape: (nz, nx)

    x_mid = 0.5 * (x_bins[:-1] + x_bins[1:])
    z_mid = 0.5 * (z_bins[:-1] + z_bins[1:])

    # ── 6. interface detection — scan each z-row, find left & right boundary ──
    # density shape: (nz, nx)  →  row = density[iz, :]
    dx = float(x_bins[1] - x_bins[0])
    left_x, left_z, right_x, right_z = [], [], [], []

    for iz, z_val in enumerate(z_mid):
        if z_val < Z_FIT_MIN:
            continue
        row = density[iz, :]
        above = np.where(row >= rho_iface)[0]
        if len(above) == 0:
            continue

        li = int(above[0])
        ri = int(above[-1])

        # subpixel left boundary
        if li > 0:
            denom = row[li] - row[li-1]
            lx = x_mid[li-1] + ((rho_iface - row[li-1]) / denom * dx
                                 if abs(denom) > 1e-30 else 0.0)
        else:
            lx = x_mid[li]

        # subpixel right boundary
        if ri + 1 < len(x_mid):
            denom = row[ri] - row[ri+1]
            rx = x_mid[ri] + ((row[ri] - rho_iface) / denom * dx
                               if abs(denom) > 1e-30 else 0.0)
        else:
            rx = x_mid[ri]

        left_x.append(lx);  left_z.append(z_val)
        right_x.append(rx); right_z.append(z_val)

    # IQR outlier removal on x-positions
    def _iqr_clean(xs, zs):
        xs, zs = np.array(xs), np.array(zs)
        if len(xs) < 4:
            return xs, zs
        q1, q3 = np.percentile(xs, [25, 75])
        iqr = q3 - q1
        ok = (xs >= q1 - 1.5*iqr) & (xs <= q3 + 1.5*iqr)
        return xs[ok], zs[ok]

    left_x,  left_z  = _iqr_clean(left_x,  left_z)
    right_x, right_z = _iqr_clean(right_x, right_z)

    print(f"  Interface points : left={len(left_x)}, right={len(right_x)}")

    if len(left_x) < MIN_IFACE_PTS or len(right_x) < MIN_IFACE_PTS:
        print(f"  [ERROR] Too few interface points (need ≥ {MIN_IFACE_PTS} per side).")
        print("  Try reducing Z_FIT_MIN or increasing SLICE_HALF_WIDTH.")
        return

    # ── 7. circle fits ────────────────────────────────────────────────────────
    all_x = np.concatenate([left_x, right_x])
    all_z = np.concatenate([left_z, right_z])

    try:
        _, cz_g, R_g, rms_g = _fit_circle(all_x, all_z)
        theta_g = _contact_angle(cz_g, R_g)
        print(f"  Global fit : cz={cz_g:.3f}  R={R_g:.3f}  θ={theta_g:.2f}°  RMS={rms_g:.4f}")
    except ValueError as e:
        print(f"  [ERROR] Global fit failed: {e}")
        return

    try:
        _, cz_L, R_L, _ = _fit_circle(left_x, left_z)
        theta_L = _contact_angle(cz_L, R_L)
    except ValueError:
        theta_L = theta_g

    try:
        _, cz_R, R_R, _ = _fit_circle(right_x, right_z)
        theta_R = _contact_angle(cz_R, R_R)
    except ValueError:
        theta_R = theta_g

    theta_avg = (theta_L + theta_R) / 2.0

    # ── 8. result ─────────────────────────────────────────────────────────────
    print(f"\n{'─'*60}")
    print(f"  θ_L   = {theta_L:.2f}°   (left  contact angle)")
    print(f"  θ_R   = {theta_R:.2f}°   (right contact angle)")
    print(f"  θ_avg = {theta_avg:.2f}°   (average)")
    print(f"{'═'*60}\n")


run()
