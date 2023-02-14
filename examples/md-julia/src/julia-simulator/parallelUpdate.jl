# include("SimulatorModules.jl")
# using .Simulator, .Simulator.Properties, .Simulator.Options, .Simulator.Particles, .Simulator.AutoPasM, .Simulator.Iterators

include("SimulatorModules.jl")
using .Simulator, .Simulator.AutoPasInterface, .Simulator.Properties, .Simulator.Iterators, .Simulator.Options, .Simulator.Particles
using Base.Threads: @spawn, @threads
# create/get InputParameters
using Profile, PProf, ProfileView

grid = CubeGridInput()

grid.particlesPerDimension = [300000, 1, 1]
grid.particleSpacing = 1.1225
grid.bottomLeftCorner = [1.5, 1.5, 1.5]
grid.velocity = [0.0, 0.0, 0.0]
grid.particleType = 0
grid.particleEpsilon = 5
grid.particleSigma = 1.0
grid.particleMass = 1.0
grid.factorBrownianMotion = 0.1

inputParameters = InputParameters()

inputParameters.container = [Options.linkedCells] # vector of contianer options -> parsing needed 1
inputParameters.verletRebuildFrequency = 3
inputParameters.verletSkinRadiusPerTimestep = 1.3
inputParameters.verletClusterSize = 12
inputParameters.selectorStrategy = Options.SelectorStrategyOption(Options.fastestAbs)
inputParameters.dataLayout = [Options.aos]
inputParameters.traversal = [Options.lc_c01]
inputParameters.tuningStrategy = Options.TuningStrategyOption(Options.fullSearch)
inputParameters.mpiStrategy = Options.MPIStrategyOption(Options.noMPI)
inputParameters.tuningInterval = 12
inputParameters.tuningSamples = 12
inputParameters.tuningMaxEvidence = 3
inputParameters.functor = "strategy"  # functor option e.g Lennard-Jones (12-6) 7
inputParameters.newton3 = [Options.disabled, Options.enabled]
inputParameters.cutoff = 10
inputParameters.boxMin = [-10000000.0, -10000000.0, -10000000.0]
inputParameters.boxMax = [10000000.0, 10000000.0, 10000000.0]
inputParameters.cellSize = [1] # check which values can be used
inputParameters.deltaT = 0.0002
inputParameters.iterations = 10 # 10000000
inputParameters.globalForce = [0.0, 0.0, 1.3]
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

# start simulation: calculate new positions, forces and velocities
function updating(autoPasContainer, particlePropertiesLibrary, inputParameters)
    for iteration = 1:inputParameters.iterations
        updatePositions(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary, inputParameters.globalForce)
    end
end
for i in 1:5
    # @time updating(autoPasContainer, particlePropertiesLibrary, inputParameters)
end

function doThePrint(particle)
    println(Particles.toString(particle))
end

function particlePrinter(autoPasContainer)
    iter = AutoPasInterface.begin(autoPasContainer, Options.IteratorBehavior(Options.ownedOrHalo))
    println("in particlePrinter")
    @spawn while isValid(iter)
        println("in while loop")
        doThePrint(Simulator.Iterators.:*(iter))
        Simulator.Iterators.:++(iter)
    end
end
function updateParticle(particleRef)
    particle = Simulator.Iterators.:*(particleRef)
    println("particle id: ", Particles.getID(particle))
    println("tid: ", Threads.threadid())
end

function parallelForLoop(autoPasContainer)
    println("starting parallel")
    nP = 10
    iter = AutoPasInterface.begin(autoPasContainer, Options.IteratorBehavior(Options.ownedOrHalo))
    @threads for i in 1:10
        updateParticle(iter)
        Simulator.Iterators.:++(iter)
    end
end

function update2(autoPasContainer, iterator, startPos, numberParticles)
    # println("tid: " * string(Threads.threadid()) * " startPos: " * string(startPos))
    index = 1
    while (index < startPos) && (isValid(iterator))
        # println(index)
        index += 1
       Simulator.Iterators.:++(iterator) 
    end
    index = 1
    while (index <= numberParticles) && (AIterators.isValid(iterator))
        # println("particleID: ", Particles.getID(AIterators.:*(iterator)), " pid: ", Threads.threadid())
        index += 1
        particle = Simulator.Iterators.:*(iterator)
        Particles.addPosition(particle, [1.0, 2.0, 3.0])
        Simulator.Iterators.:++(iterator)
    end
end

function parallelUpdatePosition(autoPasContainer, deltaT, particlePropertiesLibrary, globalForce, iter, startPos, numberParticles)

    index = 1
    while (index < startPos) && (isValid(iter))
        index += 1
        Simulator.Iterators.:++(iter) 
    end

    index = 1
    while isValid(iter) && (index < numberParticles)
        particle = Simulator.Iterators.:*(iter)
        velocity = Particles.getVelocity(particle)
        force = Particles.getForce(particle)
        Particles.setOldF(particle, force)
        Particles.setForce(particle, globalForce)
        v = velocity * deltaT
        f = force * (deltaT * deltaT / (2*Properties.getMass(particlePropertiesLibrary, 0)))
        # f = force * (deltaT * deltaT / (2*Properties.getMass(particlePropertiesLibrary, Particles.getTypeId(particle))))
        # TODO: change back
        # if Particles.getID(particle) > 3
        #     Particles.addPos(particle, [100.0, 100.0, 100.0])
        # end
        Particles.addPosition(particle, v + f)
        Simulator.Iterators.:++(iter)
        index += 1
    end
end

function parallelUpdateVelocity(autoPasContainer, deltaT, particlePropertiesLibrary, iter, startPos, numberParticles)

    index = 1
    while (index < startPos) && (isValid(iter))
        index += 1
        Simulator.Iterators.:++(iter) 
    end

    index = 1
    while isValid(iter) && (index < numberParticles)

        particle = Simulator.Iterators.:*(iter)
        force = Particles.getForce(particle)
        oldForce = Particles.getOldF(particle)
        newV = (oldForce + force) * (deltaT / 2 * Properties.getMass(particlePropertiesLibrary, 0))
        # newV = (oldForce + force) * (deltaT / 2 * Properties.getMass(particlePropertiesLibrary, Particles.getTypeId(particle)))
        Particles.addVelocity(particle, newV)
        Simulator.Iterators.:++(iter)
        index += 1

    end

end

function parallelUpdate(autoPasContainer, numberThreads)
    # numberThreads = 2
    numberParticles = 1000000
    particleTmp = numberParticles / numberThreads
    startPos = 1
    iteration = 0
    while iteration < inputParameters.iterations

        iteration += 1
        @threads for i in 1:numberThreads
            iter = AutoPasInterface.begin(autoPasContainer, Options.IteratorBehavior(Options.ownedOrHalo))
            parallelUpdatePosition(autoPasContainer, 0.001, particlePropertiesLibrary, [0.0, 0.0, 1.3], iter, startPos + (i-1)*particleTmp, particleTmp)
        end

        updateForces(autoPasContainer, [0.0, 0.0, 1.3], particlePropertiesLibrary)

        @threads for i in 1:numberThreads
            iter = AutoPasInterface.begin(autoPasContainer, Options.IteratorBehavior(Options.ownedOrHalo))
            parallelUpdateVelocity(autoPasContainer, 0.001, particlePropertiesLibrary, iter, startPos + (i-1) * particleTmp, particleTmp)
        end
    end
    totalParticles = AutoPasInterface.getNumberOfParticles(autoPasContainer, Options.IteratorBehavior(Options.ownedOrHalo))
    println("total particles: " * string(totalParticles))
end

function parallelUpdatePositionBase(autoPasContainer, numberThreads, numberParticles)
    # numberParticles = 1000000
    particleTmp = numberParticles / numberThreads
    startPos = 1
    iteration = 0
    while iteration < inputParameters.iterations

        iteration += 1
        @threads for i in 1:numberThreads
            iter = AutoPasInterface.begin(autoPasContainer, Options.IteratorBehavior(Options.ownedOrHalo))
            parallelUpdatePosition(autoPasContainer, 0.001, particlePropertiesLibrary, [0.0, 0.0, 1.3], iter, startPos + (i-1)*particleTmp, particleTmp)
        end
    end
end

obj = inputParameters.objects[1].particlesPerDimension
nPar = obj[1] * obj[2] * obj[3]
println("cla part number: " * string(nPar))

println("################")
println("option 1")
for i in 1:5
    # @time parallelUpdate(autoPasContainer, 2)
    @time parallelUpdatePositionBase(autoPasContainer, Threads.nthreads(), nPar)
end

autoPasContainer2, particlePropertiesLibrary2 = parseInput(inputParameters)

function parallelUpdatePosition2(aPC, deltaT, ppL, globalForce, iter, nt)

    while isValid(iter)
        particle = Simulator.Iterators.:*(iter)
        velocity = Particles.getVelocity(particle)
        force = Particles.getForce(particle)
        setOldF(particle, force)
        setForce(particle, globalForce)
        velocity = velocity .* deltaT
        force = force .* (deltaT * deltaT / (2*Properties.getMass(particlePropertiesLibrary2, 0)))
        addPosition(particle, velocity .+ force )

        for i in 1:nt
            Simulator.Iterators.:++(iter)
        end
    end
end

function parallelUpdate2(aPC, numberThreads, pPL)
    iteration = 1

    while iteration < inputParameters.iterations
        @threads for i in 1:numberThreads
            iter = AutoPasInterface.begin(aPC, IteratorBehavior(ownedOrHalo))
            for j in 1:(i-1)
                Simulator.Iterators.:++(iter)
            end
            parallelUpdatePosition2(aPC, 0.001, pPL, [0.0, 0.0, 1.3], iter, numberThreads)
        end
        iteration += 1
    end
end


println("########################")
println("option 2")
for i in 1:5
    @time parallelUpdate2(autoPasContainer2, Threads.nthreads(), particlePropertiesLibrary2)
end

function moreCall(i)

    for j in 1:i
        parallelUpdate2(autoPasContainer2, Threads.nthreads(), particlePropertiesLibrary2)
    end
end

# ProfileView.set_graphtype!(:icicle)

# ProfileView.@profview moreCall(30)

# Profile.print()
#=
for i in 1:2
    @profile parallelUpdate2(autoPasContainer2, Threads.nthreads(), particlePropertiesLibrary2)
    Profile.print()
    # alternative 2
    # ProfileView.@profview parallelUpdate2(autoPasContainer2, Threads.nthreads(), particlePropertiesLibrary2)
end
=#

# @trace parallelUpdate2(autoPasContainer2, Threads.nthreads(), particlePropertiesLibrary2)

println("option 3")

function parallelUpdate3(aPC, numberThreads, pPL)
    iteration = 1

    while iteration < inputParameters.iterations
        @sync for i in 1:numberThreads
            iter = AutoPasInterface.begin(aPC, IteratorBehavior(ownedOrHalo))
            for j in 1:(i-1)
                Simulator.Iterators.:++(iter)
            end
            @spawn parallelUpdatePosition2(aPC, 0.001, pPL, [0.0, 0.0, 1.3], iter, numberThreads)
        end
        iteration += 1
    end
end

function spawnTest(i)
    for i in 1:i
        parallelUpdate3(autoPasContainer2, Threads.nthreads(), particlePropertiesLibrary2)
    end
end

@profile spawnTest(1)

@profile spawnTest(30)
Profile.print()