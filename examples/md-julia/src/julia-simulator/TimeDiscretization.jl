"""
calculate the new position for all particles in the particle container
@param 
"""

# using ..Simulator.AutoPasInterface, ..Simulator.Iterators

function updatePositions(autoPasContainer, deltaT, particlePropertiesLibrary, globalForce)

    iter = AutoPasInterface.begin(autoPasContainer, Options.IteratorBehavior(Options.ownedOrHalo))
    while isValid(iter)
        particle = Simulator.Iterators.:*(iter)
        velocity = getVelocity(particle)
        force = getForce(particle)
        setOldF(particle, force)
        setForce(particle, globalForce)
        velocity .*=  deltaT
        force .*= (deltaT * deltaT / (2*getMass(particlePropertiesLibrary, Particles.getTypeId(particle))))
        force .+= velocity
        addPosition(particle, force)
        Simulator.Iterators.:++(iter)
    end

end

function updatePo(autoPasContainer, deltaT, particlePropertiesLibrary, globalForce)
    iter = AutoPasInterface.begin(autoPasContainer, Options.IteratorBehavior(Options.ownedOrHalo))
    while isValid(iter)
        particle = Simulator.Iterators.:*(iter)
        velocity = [getV(particle, i) for i in 0:2]
        force = [getF(particle, i) for i in 0:2]
        setOldF(particle, force)
        setForce(particle, globalForce)
        velocity .*=  deltaT
        force .*= (deltaT * deltaT / (2*getMass(particlePropertiesLibrary, Particles.getTypeId(particle))))
        force .+= velocity
        addPosition(particle, force)
        Simulator.Iterators.:++(iter)
    end

end

function updatePositionParallelSub(autoPasContainer, deltaT, particlePropertiesLibrary, globalForce, nthreads, iter)
    # println("in position parallel sub")
    while isValid(iter)
        # println("inside the loop")
        particle = Simulator.Iterators.:*(iter)
        velocity = [getV(particle, i) for i in 0:2]
        force = [getF(particle, i) for i in 0:2]
        setOldF(particle, force)
        setForce(particle, globalForce)
        velocity .*=  deltaT
        force .*= (deltaT * deltaT / (2*getMass(particlePropertiesLibrary, Particles.getTypeId(particle))))
        force .+= velocity
        addPosition(particle, force)
        # println(toString(particle))
        for i in 1:nthreads
            Simulator.Iterators.:++(iter)
        end
    end
end

function updatePositionParallel(autoPasContainer, deltaT, particlePropertiesLibrary, globalForce)
    nthreads = Threads.nthreads()
    Threads.@threads for i in 1:nthreads
        # println("pos i: ", i)
        iter = AutoPasInterface.begin(autoPasContainer, IteratorBehavior(ownedOrHalo))
        for j in 1:(i-1)
            Simulator.Iterators.:++(iter)
        end
        updatePositionParallelSub(autoPasContainer, deltaT, particlePropertiesLibrary, globalForce, nthreads, iter)
    end
end

function updatePositionJM(pV, deltaT, ppL, gloablForce)

    for particle in pV
        velocity = particle.vel
        force = particle.f

        particle.f = gloablForce
        particle.of = force

        velocity *= deltaT
        force *= (deltaT * deltaT / (2*getMass(ppL, particle.tId)))
        _pos = velocity + force
        particle.pos += _pos
    end

end

function convertAndSo(autoPasContainer, pV)

    iter = AutoPasInterface.begin(autoPasContainer, Options.IteratorBehavior(Options.ownedOrHalo))

    while isValid(iter)
        particle = Simulator.Iterators.:*(iter)

        id_ = getID(particle)

        for p in pV
            if p.mId == id_
                p.vel = getVelocity(particle)
                p.f = getForce(particle)
                break
            end

        end
        
        Simulator.Iterators.:++(iter)
    end
end

function updateVelocities(autoPasContainer, deltaT, particlePropertiesLibrary)

    iter = AutoPasInterface.begin(autoPasContainer, Options.IteratorBehavior(Options.ownedOrHalo))

    while isValid(iter)
        particle = Simulator.Iterators.:*(iter)
        force = [getF(particle, i) for i in 0:2]
        oldForce = [getOldF(particle, i) for i in 0:2]
        force .+= oldForce
        force .*= (deltaT / 2 * getMass(particlePropertiesLibrary, Particles.getTypeId(particle)))
        addVelocity(particle, force)
        Simulator.Iterators.:++(iter)
    end

end

function updateVelocitiesParallelSub(autoPascontainer, deltaT, particlePropertiesLibrary, nthreads, iter)
    # println("thread id: ", Threads.threadid())
    while isValid(iter)
        particle = Simulator.Iterators.:*(iter)
        force = [getF(particle, i) for i in 0:2]
        oldForce = [getOldF(particle, i) for i in 0:2]
        force .+= oldForce
        force .*= (deltaT / 2 * getMass(particlePropertiesLibrary, Particles.getTypeId(particle)))
        addVelocity(particle, force)
        for i in 1:nthreads
            Simulator.Iterators.:++(iter)
        end
    end 
end

function updateVelocitiesParallel(autoPasContainer, deltaT, particlePropertiesLibrary)
    nthreads = Threads.nthreads()
    Threads.@threads for i in 1:nthreads
        # println("velo i: ", i)
        iter = AutoPasInterface.begin(autoPasContainer, IteratorBehavior(ownedOrHalo))
        for j in 1:(i-1)
            Simulator.Iterators.:++(iter)
        end
        updateVelocitiesParallelSub(autoPasContainer, deltaT, particlePropertiesLibrary, nthreads, iter)
    end
end

function handleBoundaries(autoPasContainer, minBorder, maxBorder)

    iter = AutoPasInterface.begin(autoPasContainer, Options.IteratorBehavior(Options.owned))
    
    while isValid(iter)
        particle = Simulator.Iterators.:*(iter)
        pos = Particles.getPosition(particle)
        for dIndex = 1:3
            # if pos[dIndex] < minBorder[dIndex] || pos[dIndex] > maxBorder
            if pos[dIndex] > 2
                println("delete: " * Particles.toString(particle))
                deleteParticle(autoPasContainer, particle)
                
                break
            end
        end
        Simulator.Iterators.:++(iter)
    end
end

function handlePeriodic(autoPasContainer, minBorder, maxBorder)

    iter = AutoPasInterface.begin(autoPasContainer, Options.IteratorBehavior(Options.ownedOrHalo))
    
    while isValid(iter)
        particle = Simulator.Iterators.:*(iter)
        pos = getPosition(particle)

        for dIndex = 1:3
            if pos[dIndex] < minBorder[dIndex]
                pos[dIndex] = maxBorder[dIndex]
            elseif pos[dIndex] > maxBorder[dIndex]
                pos[dIndex] = minBorder[dIndex]
            end
            setPosition(particle, pos)
        end
        Simulator.Iterators.:++(iter)
    end

end

function executeUpdates()

    particle = Particles.MoleculeJ{Float64}()

    Particles.setPosition(particle, [1.0, 2.0, 3.0])
    Particles.setVelocity(particle, [1.0, 1.0, 1.0])

    aP = AutoPas{Particles.MoleculeJ{Float64}}()
    setBoxMin(aP, [0.0, 0.0, 0.0])
    setBoxMax(aP, [7.0, 7.0, 7.0])
    init(aP)

    addParticle(aP, particle)
    iO = Options.IteratorBehavior(Options.owned)
    iTer = AutoPasInterface.begin(aP, iO)
    println("##############")
    println("init configuration")
    while isValid(iTer)
        particle = Simulator.Iterators.:*(iTer)
        println(Particles.toString(particle))
        Simulator.Iterators.:++(iTer)
    end
    println("##############")
    println("init configuration")
    ppL = Properties.ParticlePropertiesLibrary{Float64, Int32}(1.5)
    Properties.addType(ppL, 0, 1.3, 1.5, 5.0)
    println("##############")
    println("done ppL")
    updatePositions(aP, 0.0153, ppL, [0.0, 0.0, 1.3])
    println("##############")
    println("done updateing")

    updateVelocities(aP, 0.01, ppL)

    iTer = AutoPasInterface.begin(aP, iO)
    println("##############")
    println("end configuration")
    while isValid(iTer)
        particle = Simulator.Iterators.:*(iTer)
        println(Particles.toString(particle))
        Simulator.Iterators.:++(iTer)
    end
end

# executeUpdates()