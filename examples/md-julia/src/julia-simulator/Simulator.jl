function startSimulation(autoPasContainer, particlePropertiesLibrary, inputParameters)
    iter = AutoPasM.begin(autoPasContainer, Options.IteratorBehavior(Options.owned))
    while AIterators.isValid(iter)
        println(Particles.toString(AIterators.:*(iter)))
        AIterators.:++(iter)
    end

    for iteration = 1 : inputParameters.iterations
        println(iteration)
        updatePositions(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary, inputParameters.globalForce)
        
        # handle boundary particles
        # handleBoundaries(autoPasContainer, inputParameters.boxMin, inputParameters.boxMax)
        # for outflow maybe only use updateContainer which has same effekt
        updateContainer(autoPasContainer)
        
        ff.updateForces(autoPasContainer, inputParameters.globalForce, particlePropertiesLibrary)

        updateVelocities(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary)
        # println("after updateVelocities")
    end
    println()
    println("print remaining particles")
    iter = AutoPasM.begin(autoPasContainer, Options.IteratorBehavior(Options.owned))
    while AIterators.isValid(iter)
        println(Particles.toString(AIterators.:*(iter)))
        AIterators.:++(iter)
    end
    println("done")
end