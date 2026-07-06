#!/usr/bin/env python3
# AutoPas/examples/md-flexible/input/surface_generator_lj_atm.py
#
# Generates matched YAML pairs for the ATM wetting study:
#   baseline  — pure LJ (pairwise functor only)
#   treatment — LJ + ATM (identical LJ part + triwise functor-3b)
#
# The LJ physics is byte-for-byte identical in both files of a pair. Any
# contact-angle difference between baseline and treatment isolates the ATM effect.
#
# Wall ATM mode controls whether wall particles enter ATM triplets:
#   uniform — NU_WALL = NU_FLUID = 0.073  (ATM acts on fluid-wall triplets; primary case)
#   off     — NU_WALL = 0                 (cbrt mixing zeroes wall-touching triplets; control case)
#
# Usage:
#   python3 surface_generator_lj_atm.py <surface|all> [uniform|off]

import sys

VALID_SURFACES = ["smooth", "boss", "pit", "grid", "dual_boss"]

BOX_SIZE_X, BOX_SIZE_Y = 160, 160
WALL_SIZE = 140
WALL_X0 = (BOX_SIZE_X - WALL_SIZE) // 2
WALL_Y0 = (BOX_SIZE_Y - WALL_SIZE) // 2
WALL_CENTER_X = WALL_X0 + WALL_SIZE // 2
WALL_CENTER_Y = WALL_Y0 + WALL_SIZE // 2

ITERATIONS = 1000000
DELTA_T = 0.0005
GRAVITY_Z = -0.0001
TEMPERATURE = 0.3
VTK_WRITE_FREQUENCY = 1000

EPSILON_FF = 1.0
EPSILON_WW = 1.0
# Fluid-wall wetting: effective_epsilon_fw = ZETA * EPSILON_FF.
# Must match FLUID_WALL_EPSILON_SCALE in SimulationParticleTypes.h.
ZETA = 0.6
RADIUS = 17
CUTOFF = 3

# Standard ATM nu for Argon-like fluids; fixed physical constant, not a tuning knob.
NU_FLUID = 0.073

DROPLET_CENTER_Z = {
    "smooth": 24,
    "boss": 34,
    "grid": 36,
    "pit": 30,
    "dual_boss": 38,
}


def block(block_id, x, y, z, sx, sy, sz, type_id):
    return f"""    {block_id}:
      particle-spacing           :  1.0
      bottomLeftCorner           :  [{x}, {y}, {z}]
      box-length                 :  [{sx}, {sy}, {sz}]
      velocity                   :  [0, 0, 0]
      particle-type-id           :  {type_id}

"""


def generate_surface_blocks(surface_type):
    blocks = ""
    block_id = 0

    x0, y0 = WALL_X0, WALL_Y0
    wall_size = WALL_SIZE
    wall_x_end = x0 + wall_size
    wall_y_end = y0 + wall_size

    period = 6
    feature_size = 3
    feature_height = 8

    base_z = 1
    base_height = 4
    z_top = base_z + base_height

    def feature_positions(start, end, size):
        last_start = end - size
        positions = list(range(start, last_start + 1, period))
        if positions[-1] != last_start:
            positions.append(last_start)
        return positions

    if surface_type in ["smooth", "boss", "grid", "dual_boss"]:
        blocks += block(block_id, x0, y0, base_z, wall_size, wall_size, base_height, 1)
        block_id += 1

    if surface_type == "smooth":
        return blocks

    if surface_type == "boss":
        for x in feature_positions(x0, wall_x_end, feature_size):
            for y in feature_positions(y0, wall_y_end, feature_size):
                blocks += block(block_id, x, y, z_top, feature_size, feature_size, feature_height, 2)
                block_id += 1

    elif surface_type == "pit":
        lower_base_height = 2
        top_layer_z = base_z + lower_base_height
        top_layer_height = 6

        blocks += block(block_id, x0, y0, base_z, wall_size, wall_size, lower_base_height, 3)
        block_id += 1

        tile = 2

        for x in range(x0, wall_x_end, tile):
            for y in range(y0, wall_y_end, tile):
                x_mod = (x - x0) % period
                y_mod = (y - y0) % period

                is_pit = x_mod < feature_size and y_mod < feature_size

                if not is_pit:
                    blocks += block(block_id, x, y, top_layer_z, tile, tile, top_layer_height, 3)
                    block_id += 1

    elif surface_type == "grid":
        tile = 2

        for x in range(x0, wall_x_end, tile):
            for y in range(y0, wall_y_end, tile):
                x_mod = (x - x0) % period
                y_mod = (y - y0) % period

                is_grid_ridge = x_mod < feature_size or y_mod < feature_size

                if is_grid_ridge:
                    blocks += block(block_id, x, y, z_top, tile, tile, feature_height, 4)
                    block_id += 1

    elif surface_type == "dual_boss":
        lower_size = 4
        lower_height = 6

        upper_size = 2
        upper_height = 4

        for x in feature_positions(x0, wall_x_end, lower_size):
            for y in feature_positions(y0, wall_y_end, lower_size):
                blocks += block(block_id, x, y, z_top, lower_size, lower_size, lower_height, 5)
                block_id += 1

                blocks += block(block_id, x + 1, y + 1, z_top + lower_height, upper_size, upper_size, upper_height, 5)
                block_id += 1

    return blocks


def _site_entry(type_id, epsilon, nu):
    s = f"  {type_id}:\n"
    s += f"    epsilon                      :  {epsilon}\n"
    s += f"    sigma                        :  1.\n"
    s += f"    mass                         :  1.\n"
    if nu is not None:
        s += f"    nu                           :  {nu}\n"
    s += "\n"
    return s


def _sites_yaml(is_treatment, nu_wall):
    fluid_nu = NU_FLUID if is_treatment else None
    wall_nu = nu_wall if is_treatment else None

    out = (
        f"  # Type 0: fluid. Fluid-wall LJ effective epsilon = ZETA*EPS"
        f" = {ZETA * EPSILON_FF} (overridden in C++).\n"
    )
    out += _site_entry(0, EPSILON_FF, fluid_nu)
    out += "  # Types 1-5: wall (frozen by integrator). EPSILON_WW does not control fluid-wall attraction.\n"
    if is_treatment:
        out += f"  # ATM fluid-fluid-wall triplets: nu_mixed = cbrt({NU_FLUID}^2 * {nu_wall}).\n"
    for t in range(1, 6):
        out += _site_entry(t, EPSILON_WW, wall_nu)
    return out


def _make_run_name(surface, wall_mode, is_treatment):
    tag = "lj_atm" if is_treatment else "lj"
    name = (
        f"{tag}_{surface}"
        f"_r{RADIUS}"
        f"_zeta{ZETA}"
        f"_eps{EPSILON_FF}"
        f"_temp{int(TEMPERATURE * 10):02d}"
        f"_g1e4{GRAVITY_Z}"
    )
    if is_treatment:
        name += f"_nu{NU_FLUID}_{wall_mode}"
    return name


def _make_filename(surface, wall_mode, is_treatment):
    tag = "lj_atm" if is_treatment else "lj"
    name = (
        f"sessileDrop_{tag}_{surface}"
        f"_zeta{ZETA}_eps{EPSILON_FF}"
        f"_temp{int(TEMPERATURE * 10):02d}_g{GRAVITY_Z}"
    )
    if is_treatment:
        name += f"_nu{NU_FLUID}_{wall_mode}"
    return name + ".yaml"


def generate_yaml(surface, wall_mode, is_treatment):
    surface_blocks = generate_surface_blocks(surface)
    nu_wall = NU_FLUID if wall_mode == "uniform" else 0.0
    run_name = _make_run_name(surface, wall_mode, is_treatment)
    sites = _sites_yaml(is_treatment, nu_wall)

    if is_treatment:
        header = (
            f"# LJ + ATM treatment — {surface} surface, wall mode: {wall_mode}.\n"
            f"# ATM is a pure correction layered on top of the identical LJ baseline. LJ still does\n"
            f"# all cohesion (ATM has no attractive well); ATM shifts the free-energy balance at\n"
            f"# interfaces, which we isolate by comparing contact angles against the baseline.\n"
            f"# NU_WALL = {nu_wall}: "
            + (
                "ATM acts on fluid-wall triplets (primary physical case).\n"
                if wall_mode == "uniform"
                else "wall zeroed from all ATM triplets via cbrt mixing (solid-liquid stays pure LJ).\n"
            )
        )
        atm_section = (
            "\n"
            "functor-3b                       :  axilrod-teller-muto\n"
            "traversal-3b                     :  [lc_c01]\n"
            "data-layout-3b                   :  [AoS]\n"
            "newton3-3b                       :  [disabled]\n"
        )
    else:
        header = (
            f"# LJ baseline — {surface} surface.\n"
            f"# Pure pairwise Lennard-Jones; no 3-body correction.\n"
            f"# Compare contact angle against the matched LJ+ATM treatment to isolate the ATM effect.\n"
        )
        atm_section = ""

    yaml_content = (
        header + "\n"
        f"container                        :  [LinkedCells]\n"
        f"verlet-rebuild-frequency         :  10\n"
        f"verlet-skin-radius               :  1.0\n"
        f"verlet-cluster-size              :  4\n"
        f"selector-strategy                :  Fastest-Absolute-Value\n"
        f"data-layout                      :  [AoS]\n"
        f"traversal                        :  [lc_c01]\n"
        f"tuning-strategies                :  []\n"
        f"tuning-interval                  :  2500\n"
        f"tuning-samples                   :  3\n"
        f"tuning-max-evidence              :  10\n"
        f"\n"
        f"functor                          :  Lennard-Jones\n"
        f"newton3                          :  [disabled]\n"
        f"cutoff                           :  {CUTOFF}\n"
        + atm_section
        + f"\n"
        f"box-min                          :  [0, 0, 0]\n"
        f"box-max                          :  [{BOX_SIZE_X}, {BOX_SIZE_Y}, 70]\n"
        f"cell-size                        :  [1]\n"
        f"\n"
        f"deltaT                           :  {DELTA_T}\n"
        f"iterations                       :  {ITERATIONS}\n"
        f"\n"
        f"energy-sensor                    :  rapl\n"
        f"boundary-type                    :  [reflective, reflective, reflective]\n"
        f"globalForce                      :  [0, 0, {GRAVITY_Z}]\n"
        f"\n"
        f"zeta                             :  {ZETA}\n"
        f"\n"
        f"Sites:\n"
        + sites
        + f"Objects:\n"
        f"  CubeClosestPacked:\n"
        + surface_blocks
        + f"\n"
        f"  Sphere:\n"
        f"    0:\n"
        f"      center                     :  [{WALL_CENTER_X}, {WALL_CENTER_Y}, {DROPLET_CENTER_Z[surface]}]\n"
        f"      radius                     :  {RADIUS}\n"
        f"      particle-spacing           :  1.122462048\n"
        f"      velocity                   :  [0, 0, 0]\n"
        f"      particle-type-id           :  0\n"
        f"\n"
        f"thermostat:\n"
        f"   initialTemperature            :  {TEMPERATURE}\n"
        f"   targetTemperature             :  {TEMPERATURE}\n"
        f"   deltaTemperature              :  0.1\n"
        f"   thermostatInterval            :  10\n"
        f"   addBrownianMotion             :  false\n"
        f"\n"
        f"vtk-filename                     :  {run_name}\n"
        f"vtk-write-frequency              :  {VTK_WRITE_FREQUENCY}\n"
        f"vtk-output-folder                :  {run_name}\n"
        f"no-end-config                    :  true\n"
        f"log-level                        :  info\n"
    )

    return yaml_content


def emit_pair(surface, wall_mode):
    bl_yaml = generate_yaml(surface, wall_mode, is_treatment=False)
    tr_yaml = generate_yaml(surface, wall_mode, is_treatment=True)

    bl_file = _make_filename(surface, wall_mode, is_treatment=False)
    tr_file = _make_filename(surface, wall_mode, is_treatment=True)

    with open(bl_file, "w") as f:
        f.write(bl_yaml)
    with open(tr_file, "w") as f:
        f.write(tr_yaml)

    print(f"  baseline  : {bl_file}")
    print(f"  treatment : {tr_file}")

    return bl_file, tr_file


def parse_args():
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python3 surface_generator_lj_atm.py <surface|all> [uniform|off]")
        sys.exit(1)

    surface_arg = sys.argv[1]
    wall_mode = sys.argv[2] if len(sys.argv) == 3 else "uniform"

    if wall_mode not in ("uniform", "off"):
        print(f"Unknown wall mode: {wall_mode!r}. Choose 'uniform' or 'off'.")
        sys.exit(1)

    if surface_arg == "all":
        surfaces = VALID_SURFACES
    elif surface_arg in VALID_SURFACES:
        surfaces = [surface_arg]
    else:
        print(f"Unknown surface: {surface_arg!r}. Valid: {VALID_SURFACES} or 'all'.")
        sys.exit(1)

    return surfaces, wall_mode


if __name__ == "__main__":
    surfaces, wall_mode = parse_args()
    print(f"Wall ATM mode: {wall_mode}")
    for surface in surfaces:
        print(f"Surface: {surface}")
        emit_pair(surface, wall_mode)
