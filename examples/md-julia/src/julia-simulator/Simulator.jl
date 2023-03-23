using TimerOutputs, BenchmarkTools
function updateDummy(autoPasContainer, shift)
    iter = Simulator.AutoPasInterface.begin(autoPasContainer, IteratorBehavior(ownedOrHalo))
    
    while isValid(iter)
        addPosition(Simulator.Iterators.:*(iter), shift)
        Simulator.Iterators.:++(iter)
    end
end
function startSimulation(autoPasContainer, particlePropertiesLibrary, inputParameters, domain, comm)
    #=
    iter = AutoPasM.begin(autoPasContainer, Options.IteratorBehavior(Options.owned))
    while AIterators.isValid(iter)
        println(Particles.toString(AIterators.:*(iter)))
        AIterators.:++(iter)
    end
    =#
    tmp = inputParameters.objects[1].particlesPerDimension
    maxParticles = tmp[1] * tmp[2] * tmp[3]
    tmpParticles = maxParticles
    # @inbounds for iteration = 1 : inputParameters.iterations
    for iteration = 1 : inputParameters.iterations
        # println("pos update")
        # @timeit "pos update" updatePositions(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary, inputParameters.globalForce)
        updatePositions(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary, inputParameters.globalForce)
        # updateDummy(autoPasContainer, shift)
        # deleteHaloParticlesLocalBounds(autoPasContainer, domain)
        # delete old halo Particles
        # updateContainer(autoPasContainer)

        # migrate particles
        # exchangeMigratingParticles(autoPasContainer, domain, comm)

        # apply periodic boundary condition
        # applyPeriodicBoundary(autoPasContainer, domain)

        # exchange halo particles
        # exchangeHaloParticles(autoPasContainer, domain, comm)

        # handle boundary particles
        # handleBoundaries(autoPasContainer, inputParameters.boxMin, inputParameters.boxMax)
        # for outflow maybe only use updateContainer which has same effekt
        # deleted_items = updateContainer(autoPasContainer)
        # ff.handlePeriodic(autoPasContainer, inputParameters.boxMin, inputParameters.boxMax)
        
        # @timeit "force update" updateForces(autoPasContainer, inputParameters.globalForce, particlePropertiesLibrary)
        updateForces(autoPasContainer, inputParameters.globalForce, particlePropertiesLibrary)
        
        # @timeit "vel update" updateVelocities(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary)
        updateVelocities(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary)
        # println("after updateVelocities")

        # if iteration == 2
        #     m = MoleculeJ{Float64}([0.5, 0.5, 0.5], [0.2, 0.15, 0.11], 8, 0)
        #     addParticle(autoPasContainer, m)
        # end
        
    end
    # print_timer()
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

function startSimulationEx(autoPasContainer, particlePropertiesLibrary, inputParameters, domain, comm, pV)
    #=
    println("iteration 0: ")
            iter = AutoPasInterface.begin(autoPasContainer, Options.IteratorBehavior(Options.ownedOrHalo))
            while isValid(iter)
                particle = Simulator.Iterators.:*(iter)
                if getID(particle) < 1
                    println(getVelocity(particle))
                    println(getPosition(particle))
                end
                Simulator.Iterators.:++(iter)
            end
    =#
    tmp = inputParameters.iterations / 10
    println(tmp)
    println("#threads: ", Threads.nthreads())
    index = 0
    iter = AutoPasInterface.begin(autoPasContainer, IteratorBehavior(ownedOrHalo))
    #=
    while isValid(iter) && (index < 10)
        p = Simulator.Iterators.:*(iter)
        if getID(p) < 10
            println(toString(p))
        end
        Simulator.Iterators.:++(iter)
        index += 1
    end
    =#
    for i in 1:10
        p = Simulator.Iterators.:*(iter)
        println(toString(p))
        Simulator.Iterators.:++(iter)
    end
    @inbounds for iteration = 1 : inputParameters.iterations
        if mod(iteration, tmp) == 0
            println("at iteration: ", iteration)
            println("#p: ", getNumberOfParticles(autoPasContainer, IteratorBehavior(ownedOrHalo)))
        end
    # for iteration = 1 : inputParameters.iterations
        #=
        if iteration % 10 == 0 && iteration < 100
            println("iteration ", iteration, ": ")
            iter = AutoPasInterface.begin(autoPasContainer, Options.IteratorBehavior(Options.ownedOrHalo))
            while isValid(iter)
                particle = Simulator.Iterators.:*(iter)
                if getID(particle) < 1
                    println(getVelocity(particle))
                    println(getPosition(particle))
                end
                Simulator.Iterators.:++(iter)
            end
        end

        if iteration == 2
            println("iteration ", iteration, ": ")
            iter = AutoPasInterface.begin(autoPasContainer, Options.IteratorBehavior(Options.ownedOrHalo))
            while isValid(iter)
                particle = Simulator.Iterators.:*(iter)
                if getID(particle) < 1
                    println(getVelocity(particle))
                    println(getPosition(particle))
                end
                Simulator.Iterators.:++(iter)
            end
        end
        =#
        @timeit "pos update julia" updatePositionJM(pV, inputParameters.deltaT, particlePropertiesLibrary, inputParameters.globalForce)
        
        # updatePositionJM(pV, inputParameters.deltaT, particlePropertiesLibrary, inputParameters.globalForce)
        
        # @timeit "pos update" updatePositions(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary, inputParameters.globalForce)
        # updatePositions(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary, inputParameters.globalForce)

        # @timeit "pos new" updatePo(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary, inputParameters.globalForce)
        @timeit "pos update parallel" updatePositionParallel(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary, inputParameters.globalForce)
        # updatePo(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary, inputParameters.globalForce)
        updateContainer(autoPasContainer)
        @timeit "force update" updateForces(autoPasContainer, inputParameters.globalForce, particlePropertiesLibrary)
        # updateForces(autoPasContainer, inputParameters.globalForce, particlePropertiesLibrary)
        
        # @timeit "vel update" updateVelocities(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary)
        @timeit "vel update parallel" updateVelocitiesParallel(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary)
        # updateVelocities(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary)
        
        # @timeit "convert" convertAndSo(autoPasContainer, pV)

    end
    iter = AutoPasInterface.begin(autoPasContainer, IteratorBehavior(ownedOrHalo))
    index = 0
    while isValid(iter) && (index < 10)
        p = Simulator.Iterators.:*(iter)
        println(toString(p))
        Simulator.Iterators.:++(iter)
        index += 1
    end

    println("left particles: ", getNumberOfParticles(autoPasContainer, IteratorBehavior(ownedOrHalo)))
    print_timer()
    #=
    iter = AutoPasInterface.begin(autoPasContainer, Options.IteratorBehavior(Options.ownedOrHalo))
    
    while isValid(iter)
        particle = Simulator.Iterators.:*(iter)
        if getID(particle) < 10
            println(getVelocity(particle))
            println(getPosition(particle))
        end
        Simulator.Iterators.:++(iter)
    end
    =#
end
#=

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

function simulate(autoPasContainer, particlePropertiesLibrary, inputParameters, domain, comm, shift)
    index = 1
    for iteration = 1 : inputParameters.iterations
        println("#######################################")
        # updatePositions(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary, inputParameters.globalForce)
        updateDummy(autoPasContainer, shift)

        # delte all halo particles
        # updateContainer(autoPasContainer)

        # migrate particles
        exchangeMigratingParticles(autoPasContainer, domain, comm)
        # apply periodic boundary condition
        index = applyPeriodicBoundary(autoPasContainer, domain, index)

        # exchange halo particles
        exchangeHaloParticles(autoPasContainer, domain, comm)

        updateForces(autoPasContainer, inputParameters.globalForce, particlePropertiesLibrary)
        
        updateVelocities(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary)
        # printAllParticles(autoPasContainer, comm)
    end

end

function printSimulation()
    println("in simulation")
end
=#