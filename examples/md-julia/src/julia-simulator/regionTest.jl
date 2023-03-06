include("SimulatorModules.jl")
using .Simulator.AutoPasInterface, .Simulator.Properties, .Simulator.Iterators, .Simulator.Options, .Simulator.Particles, .Simulator
using MPI
# create/get InputParameters

function printAllParticles(autoPasContainer, comm)
    ost = ""
    iter = Simulator.AutoPasInterface.begin(autoPasContainer, Options.IteratorBehavior(Options.ownedOrHalo))

    while isValid(iter)
        particle = Simulator.Iterators.:*(iter)
        ost = ost*toString(particle)*"\n"
        # println(toString(particle))
        Simulator.Iterators.:++(iter)
    end
    r = MPI.Comm_rank(comm)
    println("\n rank $r:\n $ost")
end

function init()
    grid = CubeGridInput()

    grid.particlesPerDimension = [0, 0, 0]
    grid.particleSpacing = 1.1225
    grid.bottomLeftCorner = [2.5, 1.5, 1.5]
    grid.velocity = [0.0, 0.0, 0.0]
    grid.particleType = 0
    grid.particleEpsilon = 5
    grid.particleSigma = 1.0
    grid.particleMass = 1.0
    grid.factorBrownianMotion = 0.1

    grid1 = CubeGridInput()

    grid1.particlesPerDimension = [0, 0, 0]
    grid1.particleSpacing = 1.1225
    grid1.bottomLeftCorner = [502.5, 1.5, 1.5]
    grid1.velocity = [0.0, 0.0, 0.0]
    grid1.particleType = 0
    grid1.particleEpsilon = 5
    grid1.particleSigma = 1.0
    grid1.particleMass = 1.0
    grid1.factorBrownianMotion = 0.1

    inputParameters = InputParameters()

    inputParameters.container = [linkedCells] # vector of contianer options -> parsing needed 1
    inputParameters.verletRebuildFrequency = 1
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
    inputParameters.cutoff = 1.5
    inputParameters.boxMin = [0, 0, 0]
    inputParameters.boxMax = [30.0, 30.0, 30.0]
    inputParameters.cellSize = [1] # check which values can be used
    inputParameters.deltaT = 0.0002
    inputParameters.iterations = 1
    inputParameters.globalForce = [0.0, 0.0, 0.0]
    inputParameters.periodicBoundaries = true
    inputParameters.objects = [grid, grid1]
    inputParameters.thermostat = Thermostat()
    inputParameters.logLevel = "strategy" # log level maybe string # 9
    inputParameters.noFlops = true
    inputParameters.noEndConfig = true # what does this mean?
    inputParameters.noProgressBar = true# what does this mean
    inputParameters.vtkFilename = "strategy"
    inputParameters.vtkWriteFrequency = "strategy"

    println("starting simulator")
    MPI.Init()
    comm = MPI.COMM_WORLD

    # parse input, create AutoPasContainer and ParticlePropertiesLibrary
    autoPasContainer, particlePropertiesLibrary, domain = parseInput(inputParameters, comm)
    return autoPasContainer, particlePropertiesLibrary, domain, inputParameters, comm
    
    println("ending simulation")
end

autoPasContainer, particlePropertiesLibrary, domain, inputParameters, comm = init()
println("rank: ", MPI.Comm_rank(comm))

if MPI.Comm_rank(comm) == 0
    println("in here")
    pos = [16.0, 5.0, 5.0]
    m = MoleculeJ{Float64}(pos, [1.0, 1.0, 1.0], 0, 0)
    addParticle(autoPasContainer, m)
    pos = [16.0, 5.0, 5.0]
    m1 = MoleculeJ{Float64}(pos, [1.0, 1.0, 1.0], 0, 0)
    addParticle(autoPasContainer, m1)
end

minPos = [0.0, 0.0, 0.0]
maxPos = [15.0, 15.0, 15.0]
#=
iter = regionIterator(autoPasContainer, minPos, maxPos, IteratorBehavior(ownedOrHalo))
println("region")
while isValid(iter)
    println(toString(Simulator.Iterators.:*(iter)))
    Simulator.Iterators.:++(iter)
end
println("region end")
=#
iG = Simulator.AutoPasInterface.begin(autoPasContainer, IteratorBehavior(ownedOrHalo))
println("changed pos")
while isValid(iG)
    addPosition(Simulator.Iterators.:*(iG), [-3.0, 0.0, 0.0])
    # println(toString(Simulator.Iterators.:*(iG)))
    Simulator.Iterators.:++(iG)
end
println("changed pos end")

# updateContainer(autoPasContainer)
println("updated Container")
iter = regionIterator(autoPasContainer, minPos, maxPos, IteratorBehavior(ownedOrHalo))
println("region")
while isValid(iter)
    println(toString(Simulator.Iterators.:*(iter)))
    Simulator.Iterators.:++(iter)
end
println("end region")