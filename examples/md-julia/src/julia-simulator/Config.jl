"""
function creating particles, the AutoPas container and ParticlePropertyLibrary object
"""

include("./SimulatorModules.jl")
using .Simulator.Options
using .Simulator.Particles
include("./InputConfiguration.jl")
using .InputConfiguration
include("./Generator.jl")
using .Generator
function initializeSimulator(inputConfiguration)
    #=
    for particleObject in inputConfiguration.objects
        if typeof(particleObject) == "CubeGridInput"

        elseif typeof()
    end
    =#
    x = [1,2,3]
    y = inputConfiguration
    return x, y

end

println("START")

grid = InputConfiguration.CubeGridInput()

grid.particlesPerDimension = [2, 2, 2]
grid.particleSpacing = 1.5
grid.bottomLeftCorner = [1.5, 1.5, 1.5]
grid.velocity = [1.0, 1.0, 1.0]
grid.particleType = 1
grid.particleEpsilon = 1.2
grid.particleSigma = 1.2
grid.particleMass = 1.0

particles = generateCubeGrid(grid)
pa = Particles.MoleculeJ{Float64}()
for particle in particles
    println(Particles.toString(particle))
end

inputParameters = InputConfiguration.InputParameters()

inputParameters.container = [Options.directSum] # vector of contianer options -> parsing needed 1
inputParameters.verletRebuildFrequency = 10
inputParameters.verletSkinRadiusPerTimestep = 1.3
inputParameters.verletClusterSize = 12
inputParameters.selectorStrategy = Options.SelectorStrategyOption(Options.fastestAbs)
inputParameters.dataLayout = [Options.aos]
inputParameters.traversal = [Options.lc_c01]
inputParameters.tuningStrategy = Options.TuningStrategyOption(Options.randomSearch)
inputParameters.mpiStrategy = Options.MPIStrategyOption(Options.noMPI)
inputParameters.tuningInterval = 10
inputParameters.tuningSamples = 3
inputParameters.tuningMaxEvidence = 3
inputParameters.functor = "strategy"  # functor option e.g Lennard-Jones (12-6) 7
inputParameters.newton3 = [Options.enabled]
inputParameters.cutoff = 1.3
inputParameters.boxMin = [0.0, 0.0, 0.0]
inputParameters.boxMax = [7.5, 7.5, 7.5]
inputParameters.cellSize = [10.0, 10.0, 10.0] # check which values can be used
inputParameters.deltaT = 0.001
inputParameters.iterations = 100
inputParameters.periodicBoundaries = true
inputParameters.objects = [grid]
inputParameters.thermostat = InputConfiguration.Thermostat()
inputParameters.logLevel = "strategy" # log level maybe string # 9
inputParameters.noFlops = true
inputParameters.noEndConfig = true # what does this mean?
inputParameters.noProgressBar = true# what does this mean
inputParameters.vtkFilename = "strategy"
inputParameters.vtkWriteFrequency = "strategy"

# op = InputConfiguration.parseContainerOption(inputParameters)

# println(op)
g = inputParameters.objects[1]
if typeof(g) == typeof(InputConfiguration.CubeGridInput())
    println("is cubegrid")
else
    println("is not cubegrid")
end
println(typeof(g))

a, b = initializeSimulator("hi")
println(a)
println(b)