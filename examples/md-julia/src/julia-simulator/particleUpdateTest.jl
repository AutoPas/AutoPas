using CxxWrap.StdLib
function startSimulation2(autoPasContainer, particlePropertiesLibrary, inputParameters)
    
    iter = AutoPasInterface.begin(autoPasContainer, Options.IteratorBehavior(Options.owned))
    while Simulator.Iterators.isValid(iter)
        pp = Simulator.Iterators.:*(iter)
        println("type of particle init: " * string(typeof(pp)))
        # println(Simulator.Particles.toString(Simulator.Iterators.:*(iter)))
        Simulator.Iterators.:++(iter)
    end
    pt = MoleculeJ{Float64}()
    println("type of newly created particle: " * string(typeof(pt)))
    tmp = inputParameters.objects[1].particlesPerDimension
    maxParticles = tmp[1] * tmp[2] * tmp[3]
    tmpParticles = maxParticles
    for iteration = 1 : inputParameters.iterations
        
        updatePositions(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary, inputParameters.globalForce)
        
        # handle boundary particles
        # handleBoundaries(autoPasContainer, inputParameters.boxMin, inputParameters.boxMax)
        # for outflow maybe only use updateContainer which has same effekt
        deleted_items = updateContainer(autoPasContainer)
        first_particle = deleted_items[1]

        println("type of updateContainer return: " * string(typeof(deleted_items)))
        
        println("typeof first item: " * string(typeof(first_particle)))

        pos1 = getPosition(first_particle)
        println("position of first particle: " * string(pos1))
        # ff.handlePeriodic(autoPasContainer, inputParameters.boxMin, inputParameters.boxMax)
        # updateForces(autoPasContainer, inputParameters.globalForce, particlePropertiesLibrary)
        
        # updateVelocities(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary)
        # println("after updateVelocities")
        
    end
    #= 
    println()
    np = 0
    iter = AutoPasM.begin(autoPasContainer, Options.IteratorBehavior(Options.ownedOrHalo))
    while AIterators.isValid(iter)
        np += 1
        AIterators.:++(iter)
    end
    println("remaining particles ", np)
    # println("done")
    =#
end

include("SimulatorModules.jl")
using .Simulator.AutoPasInterface, .Simulator.Properties, .Simulator.Iterators, .Simulator.Options, .Simulator.Particles, .Simulator

# create/get InputParameters

grid = CubeGridInput()

grid.particlesPerDimension = [100, 1, 1]
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
inputParameters.boxMax = [1000.0, 1000.0, 1000.0]
inputParameters.cellSize = [1] # check which values can be used
inputParameters.deltaT = 0.0002
inputParameters.iterations = 1
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
startSimulation2(autoPasContainer, particlePropertiesLibrary, inputParameters)

