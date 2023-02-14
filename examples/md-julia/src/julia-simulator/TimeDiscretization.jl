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
        v = velocity * deltaT
        # f = force * (deltaT * deltaT / (2*Properties.getMass(particlePropertiesLibrary, 0)))
        f = force * (deltaT * deltaT / (2*getMass(particlePropertiesLibrary, Particles.getTypeId(particle))))
        # TODO: change back
        # if Particles.getID(particle) > 3
        #     Particles.addPos(particle, [100.0, 100.0, 100.0])
        # end
        addPosition(particle, v + f)
        Simulator.Iterators.:++(iter)
    end

end

function updateVelocities(autoPasContainer, deltaT, particlePropertiesLibrary)

    iter = AutoPasInterface.begin(autoPasContainer, Options.IteratorBehavior(Options.ownedOrHalo))

    while isValid(iter)
        particle = Simulator.Iterators.:*(iter)
        force = Particles.getForce(particle)
        oldForce = Particles.getOldF(particle)
        # newV = (oldForce + force) * (deltaT / 2 * getMass(particlePropertiesLibrary, 0))
        newV = (oldForce + force) * (deltaT / 2 * getMass(particlePropertiesLibrary, Particles.getTypeId(particle)))
        addVelocity(particle, newV)
        Simulator.Iterators.:++(iter)
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