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
    #=
    particle = MoleculeJ{Float64}([-0.5, -0.5, -0.5], [0.0,0.0,0.0], 0, 0)
    addHaloParticle(autoPasContainer, particle)
    
    printAllParticles(autoPasContainer, comm)

    updateContainer(autoPasContainer)

    printAllParticles(autoPasContainer, comm)
    =#
    # @time startSimulation(autoPasContainer, particlePropertiesLibrary, inputParameters, domain, comm)
    println("ending simulation")
end

function executeTest()
    autoPasContainer, particlePropertiesLibrary, domain, inputParameters, comm = init()

    shift = easyBoundaryPlaneY(autoPasContainer, domain, comm)

    printAllParticles(autoPasContainer, comm)
    
    simulate(autoPasContainer, particlePropertiesLibrary, inputParameters, domain, comm, shift)

    printAllParticles(autoPasContainer, comm)
end
#=
inputParameters.cutoff = 1.5
inputParameters.boxMin = [0, 0, 0]
inputParameters.boxMax = [30.0, 30.0, 30.0]
=#

function easyHaloYLower(autoPasContainer, domain, comm)
    if MPI.Comm_rank(comm) == 0
        pos = [15.0, 1.0, 10.0]
        m = MoleculeJ{Float64}(pos, [1.0, 1.0, 1.0], 0, 0)
        addParticle(autoPasContainer, m)
        return [0.0, 0.0, 0.0] 
    end
end

# check if periodic boundary works in y direction
function easyBoundaryPlaneY(autoPasContainer, domain, comm)
    if MPI.Comm_rank(comm) == 0
        pos = [15.0, 1.0, 10.0]
        m = MoleculeJ{Float64}(pos, [1.0, 1.0, 1.0], 0, 0)
        addParticle(autoPasContainer, m)
        return [0.0, -2.0, 0.0] 
    end
end

# check if periodic boundary works in z direction
function easyBoundaryPlaneZ(autoPasContainer, domain)

end

# check if periodic boundary is not applied if localMin != globalMin
function easyBoundaryPlaneXFail(autoPasContainer, domain)

end

# check if periodic boundary is applied if localMin == globalMin
function easyBoundaryXSuccess(autoPasContainer, domain)

end

# check if boundary condition works at edge e.g. y-z edge
function easyBoundaryEdge(autoPasContainer, domain)

end

# check if edge boundary fails in case of x is not globalMin/Max
function boundaryEdgeXFail()

end

# check if periodic boundary is applied if x is globalMin/Max
function boundaryEdgeXSuccess()

end

# migrate particles to left next rank (in middle of domain)
function easyMigrateLeft(autoPasContainer, domain, comm)
    if MPI.Comm_rank(comm) == 1
        pos = [9.0, 5.0, 10.0]
        m = MoleculeJ{Float64}(pos, [1.0, 1.0, 1.0], 0, 0)
        # setOwnershipState(m, ownedS)
        addParticle(autoPasContainer, m)
        return [0.0, 0.0, 0.0] 
    end
end

# migrate particles to next right rank (in middle of domain)
function easyMigrateRight(autoPasContainer, domain, comm)
    if MPI.Comm_rank(comm) == 1
        pos = [21.0, 5.0, 10.0]
        m = MoleculeJ{Float64}(pos, [1.0, 1.0, 1.0], 0, 0)
        # setOwnershipState(m, ownedS)
        addParticle(autoPasContainer, m)
        return [0.0, 0.0, 0.0] 
    end
end

# migate particles to next rank if particle is in boundary area
function migrateBoundaryPalneLeft(autoPasContainer, domain, comm)
    if MPI.Comm_rank(comm) == 0
        pos = [1.0, 5.0, 10.0]
        m = MoleculeJ{Float64}(pos, [1.0, 1.0, 1.0], 0, 0)
        # setOwnershipState(m, ownedS)
        addParticle(autoPasContainer, m)
        return [-2.0, 0.0, 0.0] 
    end
end

# migrate particle to next rank if particle is in boundary area
function migrateBoundaryPlaneRight(autoPasContainer, domain, comm)
    if MPI.Comm_rank(comm) == 2
        pos = [29.0, 5.0, 10.0]
        m = MoleculeJ{Float64}(pos, [1.0, 1.0, 1.0], 0, 0)
        # setOwnershipState(m, ownedS)
        addParticle(autoPasContainer, m)
        return [2.0, 0.0, 0.0] 
    end
end

function migrateBoundaryEdgeLeft(autoPasContainer, domain, comm)

end

function migrateBoundaryEdgeRight(autoPasContainer, domain, comm)

end

function easyHaloLeft(autoPasContainer, domain, comm)
    if MPI.Comm_rank(comm) == 1
        pos = [11.0, 5.0, 10.0]
        m = MoleculeJ{Float64}(pos, [1.0, 1.0, 1.0], 0, 0)
        # setOwnershipState(m, ownedS)
        addParticle(autoPasContainer, m)
        return [0.0, 0.0, 0.0] 
    end
end

function easyHaloRight(autoPasContainer, domain, comm)
    if MPI.Comm_rank(comm) == 1
        pos = [19.0, 5.0, 10.0]
        m = MoleculeJ{Float64}(pos, [1.0, 1.0, 1.0], 0, 0)
        # setOwnershipState(m, ownedS)
        addParticle(autoPasContainer, m)
        return [0.0, 0.0, 0.0] 
    end
end