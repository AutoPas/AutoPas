# include("SimulatorModules.jl")
# using .Simulator, .Simulator.Properties, .Simulator.Options, .Simulator.Particles, .Simulator.AutoPasM, .Simulator.Iterators

include("SimulatorModules.jl")
using .Simulator.AutoPasInterface, .Simulator.Properties, .Simulator.Iterators, .Simulator.Options, .Simulator.Particles, .Simulator
using MPI
# create/get InputParameters

grid = CubeGridInput()

grid.particlesPerDimension = [1, 1, 1]
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
inputParameters.boxMax = [10.0, 10.0, 10.0]
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

MPI.Init()

function printParticles(autoPasContainer, domain)
    ost = []
    iter = AutoPasInterface.begin(autoPasContainer, IteratorBehavior(ownedOrHalo))
    while isValid(iter)
        particle = Simulator.Iterators.:*(iter)
        push!(ost, toString(particle))
        Simulator.Iterators.:++(iter)
    end
    rank = domain.rank
    println("$rank: $ost")
end

comm = MPI.COMM_WORLD

println("starting simulator")
# parse input, create AutoPasContainer and ParticlePropertiesLibrary
autoPasContainer, particlePropertiesLibrary, domain = parseInput(inputParameters, comm)

# autopas container: add a new particle inside the next rank
println("adding particles")
if domain.rank == 0
    pos1 = [domain.localBoxMax[1] + 1, 3, 3]
    pos2 = [domain.localBoxMin[1] + 1, 3, 3]
else
    pos1 = [domain.localBoxMin[1] - 1, 3, 3]
    pos2 = [domain.localBoxMin[1] + 1, 3, 3]
end

particle = MoleculeJ{Float64}(pos1, [1.0, 1.0, 1.0], domain.rank+10, 0)
addParticle(autoPasContainer, particle)
particle = MoleculeJ{Float64}(pos2, [1.0, 1.0, 1.0], domain.rank+20, 0)
addParticle(autoPasContainer, particle)

printParticles(autoPasContainer, domain)

println("################################")

migrateParticles(autoPasContainer, domain)

printParticles(autoPasContainer, domain)

# println("number of particles in main: " * string(getNp(autoPasContainer)))
# startSimulation(autoPasContainer, particlePropertiesLibrary, inputParameters)
# start simulation: calculate new positions, forces and velocities
# for i in 1:5
    # @time startSimulation(autoPasContainer, particlePropertiesLibrary, inputParameters)
# end

# printSimulation()
