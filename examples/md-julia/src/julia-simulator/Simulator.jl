function startSimulation(autoPasContainer, particlePropertiesLibrary, inputParameters)
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
    for iteration = 1 : inputParameters.iterations
        
        updatePositions(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary, inputParameters.globalForce)
        
        # handle boundary particles
        # handleBoundaries(autoPasContainer, inputParameters.boxMin, inputParameters.boxMax)
        # for outflow maybe only use updateContainer which has same effekt
        # deleted_items = updateContainer(autoPasContainer)
        # ff.handlePeriodic(autoPasContainer, inputParameters.boxMin, inputParameters.boxMax)
        updateForces(autoPasContainer, inputParameters.globalForce, particlePropertiesLibrary)
        
        updateVelocities(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary)
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

function printSimulation()
    println("in simulation")
end