# Contact Angle Measurement — Methodology

## Overview

`contact_angle.py` measures the **static contact angle** of a sessile droplet
from AutoPas MD simulation snapshots (VTU/PVTU format) using a
**2D central cross-section + circle fit** approach.

---

## Step-by-Step Method

### 1. Particle Filtering

Droplet particles are identified by `typeId == 0`.
If typeId data is missing, particles above `surface_z + margin` are used as fallback.

Substrate (typeId 1) and boss/pillar particles (typeId 2) are excluded.

---

### 2. Droplet Centroid Detection

The centroid `(x₀, y₀)` of the droplet is computed from the mean position
of all droplet particles. This avoids hardcoding the box center, which fails
on asymmetric or off-center droplets.

---

### 3. Central 2D Cross-Section

A thin slab is extracted:

```
|y - y₀| < slice_half_width
```

This reduces the 3D problem to a 2D `(x, z)` plane through the droplet center,
analogous to taking a photograph from the side.

---

### 4. Density Field Construction

Particle positions in the slab are binned into a 2D histogram on the `(x, z)` grid
with bin size `bin_size` (default 0.5 σ). The histogram is smoothed with a
Gaussian filter (sigma = `smooth_sigma`) to remove thermal noise.

---

### 5. Liquid–Vapor Interface Detection — Gibbs Dividing Surface

The interface is located at the density midpoint:

```
ρ_interface = (ρ_liquid + ρ_vapor) / 2
```

Saturation densities for the LJTS fluid (Vrabec et al., Mol. Phys. 2006):

| T    | ρ_liquid | ρ_vapor |
|------|----------|---------|
| 0.3  | 0.840    | 0.0004  |
| 0.6  | 0.752    | 0.0120  |
| 0.7  | 0.715    | 0.0220  |
| 0.8  | 0.668    | 0.0430  |
| 1.0  | 0.530    | 0.1500  |

Temperature is auto-detected from the filename (e.g. `temp06` → T = 0.6).

For each column `x` in the density grid, the z-position where density crosses
`ρ_interface` is found via subpixel linear interpolation.

A **wall-exclusion zone** `z > z_fit_min` (default 2.0 σ) discards the first
few layers above the substrate where density oscillations from wall layering
(Becker et al., 2014) would corrupt the fit.

**IQR-based outlier removal** discards interface points that deviate more than
1.5 × IQR from the median, removing vapor-phase noise spikes.

---

### 6. Circle Fit — Young-Laplace Approximation

A sessile droplet in mechanical equilibrium has a spherical cap profile,
as required by the Young-Laplace equation for a constant-pressure interface.

The interface points `(x_i, z_i)` are fit to a circle:

```
(x - cx)² + (z - cz)² = R²
```

**Two-stage fit:**
1. **Algebraic least-squares** (Pratt/Taubin) — fast, closed-form initial estimate
2. **Levenberg-Marquardt geometric refinement** — minimizes true geometric distance
   to the circle, giving a more accurate result for noisy data

**Split fit:** separate circles are fit to the left half (`x < cx`) and right half
(`x > cx`) of the interface to capture asymmetry on structured surfaces
(boss, dual-boss, grid, pit).

---

### 7. Contact Angle Computation

From the circle geometry:

```
cos(θ) = −z_c / R
```

where `z_c` is the circle center height above the substrate (`z = 0`)
and `R` is the circle radius.

Sign convention:

| Condition   | cos(θ) | θ        | Surface type  |
|-------------|--------|----------|---------------|
| z_c > 0     | < 0    | > 90°    | Hydrophobic   |
| z_c = 0     | = 0    | = 90°    | Hemisphere    |
| z_c < 0     | > 0    | < 90°    | Hydrophilic   |

Three angles are reported per file:
- **θ_L** — left half fit
- **θ_R** — right half fit
- **θ_avg** — average of left and right

---

### 8. Output

| File | Content |
|------|---------|
| `contact_angle_results/<name>.png` | Diagnostic plot: density map, interface points, circle fit, arcs |
| `contact_angle_results/contact_angles.csv` | θ_L, θ_R, θ_avg, R, cx, cz per file |
| `contact_angle_results/contact_angles.txt` | Human-readable summary table |

---

## Key References

- **Vrabec et al.** (2006) — LJTS saturation densities.
  *Mol. Phys.* 104, 1509–1527.
- **Becker et al.** (2014) — Wall layering exclusion zone for contact angle fitting.
- **Young, T.** (1805) — Equilibrium of liquid on a solid surface (Young equation).
- **Young-Laplace equation** — Spherical cap profile for constant surface tension droplet.
