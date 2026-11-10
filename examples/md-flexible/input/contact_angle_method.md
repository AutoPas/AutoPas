# Contact-Angle Measurement — Method and Sources

## Overview

The contact angle of the sessile droplet is measured with **Method II of Ravipati et al. (2018)**:
the angle between the liquid–vapour interface and the wall is obtained from the interface normal and
extrapolated to the contact line. A **spherical-cap (circle) fit** (Method I of Ravipati et al.;
Becker et al., 2014) is computed on the same data as an independent cross-check.

The analysis is implemented in `contact_angle.py` and runs directly on the AutoPas md-flexible
VTK output (`*_Particles_*.pvtu`).

## Procedure

1. **Wall / fluid separation.** Particles with `typeId = 0` are fluid; all other types are wall.
   The height of the top wall layer defines `z = 0` (for structured surfaces: the top of the
   pillars, i.e. the apparent contact angle is measured).
2. **Droplet identification.** The droplet is the largest cluster of fluid particles according to the
   Stillinger criterion with cluster radius `r_cl = 1.5 σ` (Stillinger, 1963; as used by
   Ingebrigtsen & Toxvaerd, 2007). Vapour particles are thereby excluded.
3. **Alignment.** Each frame is shifted so that the droplet centre of mass lies at `x = y = 0`.
4. **Coarse-grained density.** The particle positions are convolved with a Gaussian kernel
   (Willard & Chandler, 2010) on a grid of spacing `0.25 σ`.
5. **Azimuthal average.** The 3D density is averaged around the droplet axis, giving `ρ(r, z)`
   (cylindrical binning as in Becker et al., 2014).
6. **Interface.** The liquid–vapour interface is the iso-density contour `ρ = (ρ_l + ρ_v) / 2`,
   where `ρ_l` and `ρ_v` are measured in the droplet core and in the far vapour, respectively.
7. **Interface normal.** The outward normal is `n = −∇ρ / |∇ρ|`. Its vertical component
   `n_z = cos β(z)` is averaged in height bins of `0.5 σ`.
8. **Extrapolation (Method II).** Data below `z_min = 2.5 σ` are discarded because of wall-induced
   layering. `cos β(z)` is fitted linearly and extrapolated to the wall:
   `cos θ = cos β(z = 0)`. A linear `cos β(z)` corresponds to a spherical cap and serves as a
   self-check of the measurement.
9. **Cross-check (Method I).** A circle is fitted to the same interface points (z ≥ z_min); the
   contact angle follows from `cos θ = −z_c / R`, with `z_c` the centre height and `R` the radius.

## Averaging and uncertainties

- **Main value:** the density fields of the **last 20 output frames** are averaged, then a single
  fit is performed (supervisor's recommendation; Becker et al. (2014) likewise determine the angle
  from the density profile averaged over the production period).
  Default frame spacing: every 10 000 steps, i.e. the last 200 000 of 1 000 000 iterations.
- **Statistical uncertainty:** θ is computed separately for each of the 20 frames; the standard
  deviation (SD) is reported.
- **Systematic uncertainty:** θ is recomputed for coarse-graining widths `ξ ∈ {1.0, 1.3, 1.6} σ`
  and exclusion heights `z_min ∈ {2.0, 2.5, 3.0, 3.5} σ`; half of the resulting spread is reported.
  This quantifies the sensitivity within the tested ranges only.

## Parameters

| Parameter | Value | Origin |
|---|---|---|
| Cluster radius `r_cl` | 1.5 σ | Stillinger (1963); Ingebrigtsen & Toxvaerd (2007) |
| Coarse-graining width `ξ` (x, y) | 1.3 σ | Ravipati et al. (2018) |
| Coarse-graining width `ξ_z` (z) | 0.5 σ | this work (see below) |
| Grid spacing | 0.25 σ | this work |
| Interface level | (ρ_l + ρ_v) / 2 | Becker et al. (2014); Ravipati et al. (2018) |
| Excluded near-wall region `z_min` | 2.5 σ | Ravipati et al. (2018) |
| Fit range of `cos β(z)` | z_min … z_min + 0.8 (z_top − z_min) | this work |
| Averaging window | last 20 frames | supervisor |

## Modifications introduced in this work

Two changes with respect to Ravipati et al. (2018) were necessary because the Gaussian kernel is cut
by the wall:

1. **Normalized convolution.** The smoothed density is divided by the smoothed mask of the region
   accessible to the fluid, so that the empty wall region is not counted as low density.
2. **Anisotropic smoothing.** Width `ξ = 1.3 σ` parallel to the wall, `ξ_z = 0.5 σ` normal to it.
   With isotropic smoothing the interface was pulled inwards near the contact line, which biased
   the angle of wetting droplets (θ < 90°) by about +2°.

### Validation with synthetic droplets of known contact angle

| True θ | Isotropic smoothing | Final method |
|---|---|---|
| 70° | 72.0° | **70.2°** |
| 95° | 95.2° | **95.3°** |
| 110° | 109.5° | **110.0°** |

## Result: LJ, smooth surface

Setup: LJ (truncated and shifted, `r_c = 2.5 σ`), ζ = 0.6, T = 0.6, 1 000 000 iterations.

| Quantity | Value |
|---|---|
| θ (Method II) | **69.8° ± 1.5° (stat., SD) ± 1.2° (syst.)** |
| θ (Method I, circle fit) | 71.1° |
| Liquid density ρ_l | 0.845 σ⁻³ |
| Vapour density ρ_v | 0.0008 σ⁻³ |

## Figures

- **`fig_profile.png`** — main figure, in the style of Becker et al. (2014).
  Left half: MD snapshot of the last frame (slice |y| < 1.5 σ through the droplet centre;
  blue = fluid, grey = wall). Right half: time-averaged coarse-grained density, fitted circle
  (dashed white), tangent at the contact line (solid white) and contact angle θ.
- **`fig_sensitivity.png`** — θ as a function of ξ and z_min (systematic uncertainty).
- With `--full`: `cos β(z)` with linear fit, θ versus simulation step, and a particle-based
  side view.

## Suggested text for the report

> The contact angle was determined following Method II of Ravipati et al. (2018), using the Gaussian
> coarse-graining of Willard and Chandler (2010) applied to the density averaged over the last 20
> frames; a spherical-cap fit (Becker et al., 2014) served as a cross-check. Fluid within 2.5 σ of
> the wall was excluded because of wall-induced layering (Ingebrigtsen and Toxvaerd, 2007).

## References

1. S. Ravipati, B. Aymard, S. Kalliadasis, A. Galindo, "On the equilibrium contact angle of sessile
   liquid drops from molecular dynamics simulations", *J. Chem. Phys.* **148**, 164704 (2018).
   https://doi.org/10.1063/1.5021088
2. A. P. Willard, D. Chandler, "Instantaneous liquid interfaces", *J. Phys. Chem. B* **114**,
   1954–1958 (2010). https://doi.org/10.1021/jp909219k
3. S. Becker, H. M. Urbassek, M. Horsch, H. Hasse, "Contact angle of sessile drops in
   Lennard-Jones systems", *Langmuir* **30**, 13606–13614 (2014).
   https://doi.org/10.1021/la503974z
4. T. Ingebrigtsen, S. Toxvaerd, "Contact angles of Lennard-Jones liquids and droplets on planar
   surfaces", *J. Phys. Chem. C* **111**, 8518–8523 (2007). https://doi.org/10.1021/jp0676235
5. F. H. Stillinger, "Rigorous basis of the Frenkel-Band theory of association equilibrium",
   *J. Chem. Phys.* **38**, 1486 (1963). https://doi.org/10.1063/1.1776907
