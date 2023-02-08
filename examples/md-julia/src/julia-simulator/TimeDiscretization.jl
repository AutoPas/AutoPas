"""
calculate the new position for all particles in the particle container
@param 
"""

function updatePositions(autoPasContainer, deltaT, particlePropertiesLibrary, globalForce)

    iter = AutoPasM.begin(autoPasContainer, Options.IteratorBehavior(Options.owned))

    while Iterators.isValid(iter)
        particle = Iterators.:*(iter)
        velocity = Particles.getV(particle)
        force = Particles.getF(particle)
        Particles.setOldF(particle, force)
        Particles.setF(particle, globalForce)
        v = velocity * deltaT
        f = force * (deltaT * deltaT / (2*Properties.getMass(particlePropertiesLibrary, 0)))
        # f = force * (deltaT * deltaT / (2*Properties.getMass(particlePropertiesLibrary, Particles.getTypeId(particle))))
        # TODO: change back
        if Particles.getID(particle) > 3
            Particles.addPos(particle, [100.0, 100.0, 100.0])
        end
        # Particles.addPos(particle, v + f)
        Iterators.:++(iter)
    end

end

function updateVelocities(autoPasContainer, deltaT, particlePropertiesLibrary)

    iter = AutoPasM.begin(autoPasContainer, Options.IteratorBehavior(Options.owned))

    while Iterators.isValid(iter)
        particle = Iterators.:*(iter)
        force = Particles.getF(particle)
        oldForce = Particles.getOldF(particle)
        newV = (oldForce + force) * (deltaT / 2 * Properties.getMass(particlePropertiesLibrary, 0))
        # newV = (oldForce + force) * (deltaT / 2 * Properties.getMass(particlePropertiesLibrary, Particles.getTypeId(particle)))
        Particles.addV(particle, newV)
        Iterators.:++(iter)
    end

end

function handleBoundaries(autoPasContainer, minBorder, maxBorder)

    iter = AutoPasM.begin(autoPasContainer, Options.IteratorBehavior(Options.owned))
    
    while Iterators.isValid(iter)
        particle = Iterators.:*(iter)
        pos = Particles.getPos(particle)
        for dIndex = 1:3
            # if pos[dIndex] < minBorder[dIndex] || pos[dIndex] > maxBorder
            if pos[dIndex] > 2
                println("delete: " * Particles.toString(particle))
                AutoPasM.deleteParticle(autoPasContainer, particle)
                
                break
            end
        end
        Iterators.:++(iter)
    end

end

function executeUpdates()

    particle = Particles.MoleculeJ{Float64}()

    Particles.setPos(particle, [1.0, 2.0, 3.0])
    Particles.setV(particle, [1.0, 1.0, 1.0])

    aP = AutoPasM.AutoPas{Particles.MoleculeJ{Float64}}()
    AutoPasM.setBoxMin(aP, [0.0, 0.0, 0.0])
    AutoPasM.setBoxMax(aP, [7.0, 7.0, 7.0])
    AutoPasM.init(aP)

    AutoPasM.addParticle(aP, particle)
    iO = Options.IteratorBehavior(Options.owned)
    iTer = AutoPasM.begin(aP, iO)
    println("##############")
    println("init configuration")
    while Iterators.isValid(iTer)
        particle = Iterators.:*(iTer)
        println(Particles.toString(particle))
        Iterators.:++(iTer)
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

    iTer = AutoPasM.begin(aP, iO)
    println("##############")
    println("end configuration")
    while Iterators.isValid(iTer)
        particle = Iterators.:*(iTer)
        println(Particles.toString(particle))
        Iterators.:++(iTer)
    end
end

# executeUpdates()