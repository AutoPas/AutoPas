# include("SimulatorModules.jl")
# using .Simulator, .Simulator.Properties, .Simulator.Options, .Simulator.Particles, .Simulator.AutoPasM, .Simulator.Iterators

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
function sim(iterations)
    grid = CubeGridInput()

    grid.particlesPerDimension = [10, 10, 10]
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
    inputParameters.cutoff = 1.5
    inputParameters.boxMin = [0, 0, 0]
    inputParameters.boxMax = [1000.0, 1000.0, 1000.0]
    inputParameters.cellSize = [1] # check which values can be used
    inputParameters.deltaT = 0.0002
    inputParameters.iterations = iterations
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
    # particle = MoleculeJ{Float64}([-0.5, -0.5, -0.5], [0.0,0.0,0.0], 0, 0)
    # addHaloParticle(autoPasContainer, particle)
    
    # printAllParticles(autoPasContainer, comm)

    # updateContainer(autoPasContainer)

    # printAllParticles(autoPasContainer, comm)
    

    @time startSimulation(autoPasContainer, particlePropertiesLibrary, inputParameters, domain, comm)
    println("ending simulation")
end
    #=
    if MPI.Comm_rank(comm) == 0
        pos = [3.7, 3.5, 3.5]
        vel = [0.01, 0.023, 0.03551]
        # p = MoleculeJ{Float64}(pos, vel, 12, 0)
        p = MoleculeJ{Float64}(pos, vel, MPI.Comm_rank(comm), 0)
        addParticle(autoPasContainer, p)

        pos = [2.7, 3.5, 3.5]
        vel = [0.01, 0.023, 0.03551]
        # p = MoleculeJ{Float64}(pos, vel, 12, 0)
        p = MoleculeJ{Float64}(pos, vel, MPI.Comm_rank(comm) + 3, 0)
        addParticle(autoPasContainer, p)

    end
    if MPI.Comm_rank(comm) == 1
        pos = [5.15, 3.5, 3.5]
        vel = [0.01, 0.023, 0.03551]
        # p = MoleculeJ{Float64}(pos, vel, 12, 0)
        p = MoleculeJ{Float64}(pos, vel, MPI.Comm_rank(comm), 0)
        addParticle(autoPasContainer, p)
    end
    =#
    # printAllParticles(autoPasContainer, comm)

    # println("############################")
    # exchangeHaloParticles(autoPasContainer, domain, comm)
    # exchangeMigratingParticles(autoPasContainer, domain, comm)
    # println("number of particles in main: " * string(getNp(autoPasContainer)))
    # startSimulation(autoPasContainer, particlePropertiesLibrary, inputParameters)
    # updateForces(autoPasContainer, inputParameters.globalForce, particlePropertiesLibrary)

    # printAllParticles(autoPasContainer, comm)

    # println("############################")

    # deleteHaloParticlesLocalBounds(autoPasContainer, domain)

    # printAllParticles(autoPasContainer, comm)

sim(1)

sim(1000)
# sim(2)
