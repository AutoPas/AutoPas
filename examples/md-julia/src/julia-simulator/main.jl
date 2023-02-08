# include("SimulatorModules.jl")
# using .Simulator, .Simulator.Properties, .Simulator.Options, .Simulator.Particles, .Simulator.AutoPasM, .Simulator.Iterators

include("SimulatorModules.jl")
using .ff, .AutoPasM, .Properties, .Iterators, .Options, .Particles 

# create/get InputParameters

grid = CubeGridInput()

grid.particlesPerDimension = [2, 2, 2]
grid.particleSpacing = 1.5
grid.bottomLeftCorner = [1.5, 1.5, 1.5]
grid.velocity = [1.0, 1.0, 1.0]
grid.particleType = 0
grid.particleEpsilon = 1.2
grid.particleSigma = 1.2
grid.particleMass = 1.0
grid.factorBrownianMotion = 0.1

inputParameters = ff.InputParameters()

inputParameters.container = [Options.directSum] # vector of contianer options -> parsing needed 1
inputParameters.verletRebuildFrequency = 3
inputParameters.verletSkinRadiusPerTimestep = 1.3
inputParameters.verletClusterSize = 12
inputParameters.selectorStrategy = Options.SelectorStrategyOption(Options.fastestAbs)
inputParameters.dataLayout = [Options.aos]
inputParameters.traversal = [Options.lc_c01]
inputParameters.tuningStrategy = Options.TuningStrategyOption(Options.randomSearch)
inputParameters.mpiStrategy = Options.MPIStrategyOption(Options.noMPI)
inputParameters.tuningInterval = 12
inputParameters.tuningSamples = 12
inputParameters.tuningMaxEvidence = 3
inputParameters.functor = "strategy"  # functor option e.g Lennard-Jones (12-6) 7
inputParameters.newton3 = [Options.enabled]
inputParameters.cutoff = 1.3
inputParameters.boxMin = [0.0, 0.0, 0.0]
inputParameters.boxMax = [17.5, 17.5, 17.5]
inputParameters.cellSize = [10.0, 10.0, 10.0] # check which values can be used
inputParameters.deltaT = 0.001
inputParameters.iterations = 100
inputParameters.globalForce = [0.0, 0.0, 1.3]
inputParameters.periodicBoundaries = true
inputParameters.objects = [grid]
inputParameters.thermostat = ff.Thermostat()
inputParameters.logLevel = "strategy" # log level maybe string # 9
inputParameters.noFlops = true
inputParameters.noEndConfig = true # what does this mean?
inputParameters.noProgressBar = true# what does this mean
inputParameters.vtkFilename = "strategy"
inputParameters.vtkWriteFrequency = "strategy"

# parse input, create AutoPasContainer and ParticlePropertiesLibrary
autoPasContainer, particlePropertiesLibrary = parseInput(inputParameters)

# start simulation: calculate new positions, forces and velocities
startSimulation(autoPasContainer, particlePropertiesLibrary, inputParameters)
