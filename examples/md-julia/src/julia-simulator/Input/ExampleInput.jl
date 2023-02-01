include("../InputConfiguration.jl")
using ..InputConfiguration

# constructing a CubeGrid with default constructor and setting
# all attributes individually
cubeGridS = CubeGrid()

cubeGridS.particles-per-dimension = [3, 3, 3]
cubeGridS.particle-spacing = 0.5
cubeGridS.bottomLeftCorner = [0,0,0]
cubeGridS.velocity = [1.0, 1.0, 1.0]
cubeGridS.particle-type = 0
cubeGridS.particle-epsilon = 1.0
cubeGridS.particle-sigma = 1.0
cubeGridS.particle-mass = 1.0

# constructing a CubeGrid with constructor
cubeGridC = CubeGrid(
    #=particles-per-dimension=# [2,3,5],
    #=particle-spacing=# 0.5,
    #=bottomLeftCorner=# [10,10,10],
    #=velocity=# [1.0, 1.0, 1.0],
    #=particle-type=# 0,
    #=particle-epsilon=# 1.0,
    #=particle-sigma=# 1.0,
    #=particle-mass=# 1,0
)

# define all other input parameters
inputParameters = InputParameters()
inputParameters.container = []
inputParameters.verlet-rebuild-frequency = 20
inputParameters.verlet-skin-radius-per-timestep = 0.01
inputParameters.verlet-cluster-size = 3
inputParameters.selector-strategy = "Fastest-Absolute-Value"
inputParameters.data-layout = ["AoS", "SoA"]
inputParameters.traversal::Vector{String}
inputParameters.tuning-strategy = "full-Search"
inputParameters.mpi-strategy = "no-mpi"
inputParameters.tuning-interval = 5000
inputParameters.tuning-samples = 3
inputParameters.tuning-max-evidence = 10
inputParameters.functor = "Lennard-Jones (12-6)"
inputParameters.newton3 = ["disabled", "enabled"]
inputParameters.cutoff = 1.0
inputParameters.box-min = [0.0, 0.0, 0.0]
inputParameters.box-max = [7.5, 7.5, 7.5]
inputParameters.cell-size = [1]
inputParameters.deltaT = 0.001
inputParameters.iterations = 10
inputParameters.periodic-boundaries = true
inputParameters.objects = [cubeGridS, cubeGridC]
inputParameters.thermostat = Thermostat()
inputParameters.log-level = "info"
inputParameters.no-flops = false
inputParameters.no-end-config = false
inputParameters.no-progress-bar = falses
