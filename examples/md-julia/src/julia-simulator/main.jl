# include("SimulatorModules.jl")
# using .Simulator, .Simulator.Properties, .Simulator.Options, .Simulator.Particles, .Simulator.AutoPasM, .Simulator.Iterators

include("SimulatorModules.jl")
using .Simulator.AutoPasInterface, .Simulator.Properties, .Simulator.Iterators, .Simulator.Options, .Simulator.Particles, .Simulator

# create/get InputParameters

grid = CubeGridInput()

grid.particlesPerDimension = [10, 8, 10]
grid.particleSpacing = 1.1225
grid.bottomLeftCorner = [1.5, 1.5, 1.5]
grid.velocity = [0.0, 0.0, 0.0]
grid.particleType = 0
grid.particleEpsilon = 5
grid.particleSigma = 1.0
grid.particleMass = 1.0
grid.factorBrownianMotion = 0.1

inputParameters = InputParameters()

inputParameters.container = [linkedCells] # vector of contianer options -> parsing needed 1
inputParameters.verletRebuildFrequency = 3
inputParameters.verletSkinRadiusPerTimestep = 1.3
inputParameters.verletClusterSize = 12
inputParameters.selectorStrategy = SelectorStrategyOption(fastestAbs)
inputParameters.dataLayout = [aos]
inputParameters.traversal = [lc_c01]
inputParameters.tuningStrategy = TuningStrategyOption(fullSearch)
inputParameters.mpiStrategy = MPIStrategyOption(noMPI)
inputParameters.tuningInterval = 12
inputParameters.tuningSamples = 12
inputParameters.tuningMaxEvidence = 3
inputParameters.functor = "strategy"  # functor option e.g Lennard-Jones (12-6) 7
inputParameters.newton3 = [disabled, enabled]
inputParameters.cutoff = 1
inputParameters.boxMin = [0, 0, 0]
inputParameters.boxMax = [100.0, 100.0, 100.0]
inputParameters.cellSize = [1] # check which values can be used
inputParameters.deltaT = 0.0002
inputParameters.iterations = 10
inputParameters.globalForce = [0.0, 0.0, 0.0]
inputParameters.periodicBoundaries = true
inputParameters.objects = [grid]
inputParameters.thermostat = Thermostat()
inputParameters.logLevel = "strategy" # log level maybe string # 9
inputParameters.noFlops = true
inputParameters.noEndConfig = true # what does this mean?
inputParameters.noProgressBar = true# what does this mean
inputParameters.vtkFilename = "strategy"
inputParameters.vtkWriteFrequency = "strategy"

println("starting simulator")
# parse input, create AutoPasContainer and ParticlePropertiesLibrary
autoPasContainer, particlePropertiesLibrary = parseInput(inputParameters)


# println("number of particles in main: " * string(getNp(autoPasContainer)))
startSimulation(autoPasContainer, particlePropertiesLibrary, inputParameters)
# start simulation: calculate new positions, forces and velocities
for i in 1:5
    # @time startSimulation(autoPasContainer, particlePropertiesLibrary, inputParameters)
end

printSimulation()
