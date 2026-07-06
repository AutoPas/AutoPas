# AutoPas/examples/md-flexible/input/surface_generator_atm.py
# Generates sessile drop YAML configs using the Axilrod-Teller-Muto (ATM) 3-body potential.
# Wall particles (type 1-5) are frozen by TimeDiscretization.cpp — no forces, no movement.
# ATM computes fluid-fluid-fluid and fluid-fluid-wall triplets only (wall-wall skipped by functor).
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
TEMPERATURE = 0.6
VTK_WRITE_FREQUENCY = 1000
RADIUS = 17

# LJ epsilon values — registered in ParticlePropertiesLibrary but NOT used by the ATM functor.
# The ATM functor only reads nu. These are here so the YAML parser does not complain.
EPSILON_FF = 1.0  # fluid-fluid (type 0)
EPSILON_WW = 1.0  # wall types (1-5), frozen — irrelevant for dynamics

# ATM parameter nu: dimensionless in reduced LJ units.
# 0.073 is the standard value for Argon-like fluids.
# Fluid-wall wetting in ATM is controlled by NU_WALL via geometric mixing:
#   nu_mixed(fluid-fluid-wall) = cbrt(NU_FLUID^2 * NU_WALL)
# There is no ZETA-style direct override for ATM (unlike LJ).
NU_FLUID = 0.073
NU_WALL = 0.04

# ATM cutoff: all three pairwise distances in a triplet must be within this value.
# Shorter than typical LJ cutoff because of the AND condition in the functor.
CUTOFF = 2.5

DROPLET_CENTER_Z = {
    "smooth": 24,
    "boss": 34,
    "grid": 36,
    "pit": 30,
    "dual_boss": 38,
}

if len(sys.argv) != 2:
    print("Usage:")
    print("  python3 surface_generator_atm.py smooth")
    print("  python3 surface_generator_atm.py boss")
    print("  python3 surface_generator_atm.py pit")
    print("  python3 surface_generator_atm.py grid")
    print("  python3 surface_generator_atm.py dual_boss")
    sys.exit(1)

surface = sys.argv[1]
RUN_NAME = (
    f"atm_{surface}"
    f"_r{RADIUS}"
    f"_nu{NU_FLUID}"
    f"_temp{int(TEMPERATURE*10):02d}"
    f"_g{GRAVITY_Z}"
)

if surface not in VALID_SURFACES:
    print(f"Unknown surface type: {surface}")
    print(f"Valid options: {VALID_SURFACES}")
    sys.exit(1)

output_filename = f"sessileDrop_atm_{surface}_nuF{NU_FLUID}_nuW{NU_WALL}_temp{int(TEMPERATURE*10):02d}_g{GRAVITY_Z}.yaml"


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


surface_blocks = generate_surface_blocks(surface)

yaml_content = f"""# This yaml file is generated automatically.
# Surface type: {surface}
# Potential: Axilrod-Teller-Muto (3-body only, no pairwise LJ)
# Wall particles (type 1-5) are frozen by the integrator — no forces applied.

container                        :  [LinkedCells]
verlet-rebuild-frequency         :  10
verlet-skin-radius               :  0.5
verlet-cluster-size              :  4
selector-strategy                :  Fastest-Absolute-Value
tuning-strategies                :  []
tuning-interval                  :  2500
tuning-samples                   :  3
tuning-max-evidence              :  10

functor-3b                       :  axilrod-teller-muto
traversal-3b                     :  [lc_c01]
newton3-3b                       :  [disabled]
data-layout-3b                   :  [AoS]
cutoff                           :  {CUTOFF}

box-min                          :  [0, 0, 0]
box-max                          :  [{BOX_SIZE_X}, {BOX_SIZE_Y}, 70]
cell-size                        :  [1]

deltaT                           :  {DELTA_T}
iterations                       :  {ITERATIONS}

energy-sensor                    :  rapl
boundary-type                    :  [reflective, reflective, reflective]
globalForce                      :  [0, 0, {GRAVITY_Z}]

Sites:
  # Type 0: fluid particles.
  # epsilon is unused by ATM functor; wetting controlled by NU_FLUID / NU_WALL.
  0:
    epsilon                      :  {EPSILON_FF}
    sigma                        :  1.
    mass                         :  1.
    nu                           :  {NU_FLUID}

  # Types 1-5: wall particles (frozen by integrator).
  # epsilon unused by ATM. nu_mixed for fluid-fluid-wall = cbrt(NU_FLUID^2 * NU_WALL).
  1:
    epsilon                      :  {EPSILON_WW}
    sigma                        :  1.
    mass                         :  1.
    nu                           :  {NU_WALL}

  2:
    epsilon                      :  {EPSILON_WW}
    sigma                        :  1.
    mass                         :  1.
    nu                           :  {NU_WALL}

  3:
    epsilon                      :  {EPSILON_WW}
    sigma                        :  1.
    mass                         :  1.
    nu                           :  {NU_WALL}

  4:
    epsilon                      :  {EPSILON_WW}
    sigma                        :  1.
    mass                         :  1.
    nu                           :  {NU_WALL}

  5:
    epsilon                      :  {EPSILON_WW}
    sigma                        :  1.
    mass                         :  1.
    nu                           :  {NU_WALL}

Objects:
  CubeClosestPacked:
{surface_blocks}

  Sphere:
    0:
      center                     :  [{WALL_CENTER_X}, {WALL_CENTER_Y}, {DROPLET_CENTER_Z[surface]}]
      radius                     :  {RADIUS}
      particle-spacing           :  1.122462048
      velocity                   :  [0, 0, 0]
      particle-type-id           :  0

thermostat:
   initialTemperature            :  {TEMPERATURE}
   targetTemperature             :  {TEMPERATURE}
   deltaTemperature              :  0.1
   thermostatInterval            :  10
   addBrownianMotion             :  false

vtk-filename                     :  {RUN_NAME}
vtk-write-frequency              :  {VTK_WRITE_FREQUENCY}
vtk-output-folder                :  {RUN_NAME}
no-end-config                    :  true
log-level                        :  info
"""

with open(output_filename, "w") as f:
    f.write(yaml_content)

print(f"Generated: {output_filename}")
