function startSimulation(autoPasContainer, particlePropertiesLibrary, inputParameters)
    iter = AutoPasM.begin(autoPasContainer, Options.IteratorBehavior(Options.owned))
    while Iterators.isValid(iter)
        println(Particles.toString(Iterators.:*(iter)))
        Iterators.:++(iter)
    end

    for iteration = 1 : inputParameters.iterations
        # println("in loop")
        updatePositions(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary, inputParameters.globalForce)
        # println("after updatePosition")
        
        # handle boundary particles
        handleBoundaries(autoPasContainer, inputParameters.boxMin, inputParameters.boxMax)
        # for outflow maybe only use updateContainer which has same effekt

        # updateForces()

        updateVelocities(autoPasContainer, inputParameters.deltaT, particlePropertiesLibrary)
        # println("after updateVelocities")
    end
    println()
    println("print remaining particles")
    iter = AutoPasM.begin(autoPasContainer, Options.IteratorBehavior(Options.owned))
    while Iterators.isValid(iter)
        println(Particles.toString(Iterators.:*(iter)))
        Iterators.:++(iter)
    end
    println("done")
end