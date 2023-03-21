include("SimulatorModule.jl")
using .Simulator.AutoPasInterface, .Simulator.Properties, .Simulator.Iterators, .Simulator.Options, .Simulator.Particles, .Simulator
using MPI, Profile, BenchmarkTools, TimerOutputs

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
    
    # define CubeGrid
    grid = CubeGridInput()

    grid.particlesPerDimension = [100, 100, 100]
    grid.particleSpacing = 1.12
    grid.bottomLeftCorner = [20.0, 20.0, 20.0]
    grid.velocity = [0.0, 0.0, 0.0]
    grid.particleType = 0
    grid.particleEpsilon = 5
    grid.particleSigma = 1.0
    grid.particleMass = 2.0
    grid.factorBrownianMotion = 0.1

    # define another CubeGrid
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

    # define a UniformGrid
    uniGrid = UniformGrid()

    uniGrid.boxLength = [98, 98, 98]
    uniGrid.numParticles = 1000
    uniGrid.bottomLeftCorner = [1.0, 1.0, 1.0]
    uniGrid.velocity = [0.0, 0.0, 0.0]
    uniGrid.particleType = 0
    uniGrid.particleEpsilon = 5.0
    uniGrid.particleSigma = 1.0
    uniGrid.particleMass = 1.0
    uniGrid.factorBrownianMotion = 0.1

    # set all InputParameters for a simulation
    inputParameters = InputParameters()

    inputParameters.container = [linkedCells] # vector of contianer options -> parsing needed 1
    inputParameters.verletRebuildFrequency = 20
    inputParameters.verletSkinRadiusPerTimestep = 0.1
    inputParameters.verletClusterSize = 4
    inputParameters.selectorStrategy = SelectorStrategyOption(fastestAbs)
    inputParameters.dataLayout = [aos]
    inputParameters.traversal = [lc_c08]
    inputParameters.tuningStrategy = TuningStrategyOption(fullSearch)
    inputParameters.mpiStrategy = MPIStrategyOption(noMPI)
    inputParameters.tuningInterval = 2500
    inputParameters.tuningSamples = 3
    inputParameters.tuningMaxEvidence = 10
    inputParameters.functor = "strategy"  # functor option e.g Lennard-Jones (12-6) 7
    inputParameters.newton3 = [disabled, enabled]
    inputParameters.cutoff = 3
    inputParameters.boxMin = [0, 0, 0]
    inputParameters.boxMax = [1500.0, 150.0, 150.0]
    inputParameters.cellSize = [1] # check which values can be used
    inputParameters.deltaT = 0.0005
    inputParameters.iterations = iterations
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
    MPI.Init()
    comm = MPI.COMM_WORLD
    reset_timer!()
    # parse input, create AutoPasContainer and ParticlePropertiesLibrary
    autoPasContainer, particlePropertiesLibrary, domain = parseInput(inputParameters, comm)

    pV, mid = generateCubeGridJulia(inputParameters.objects[1], 0) 
    # pV, mid = generateUniformGridJulia(inputParameters.objects[1], 0) 
    
    # @profile startSimulation(autoPasContainer, particlePropertiesLibrary, inputParameters, domain, comm)
    # startSimulation(autoPasContainer, particlePropertiesLibrary, inputParameters, domain, comm)

    # @profile startSimulationEx(autoPasContainer, particlePropertiesLibrary, inputParameters, domain, comm, pV)
    startSimulationEx(autoPasContainer, particlePropertiesLibrary, inputParameters, domain, comm, pV)
    
    # Profile.print()
    println("ending simulation")
end

function start()
    sim(1)
    for i in 1:1
        @time sim(3)
    end
end

start()
