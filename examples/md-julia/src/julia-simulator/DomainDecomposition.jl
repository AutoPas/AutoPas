function migrateParticles(autoPasContainer, domain)
# TODO: are halo particles set automatically? I guess they are not right

# 1. check halo particles
# 2. check particles out of boxMax
# 3. check particles out of local domain
    mp = []
    iter = AutoPasInterface.begin(autoPasContainer, IteratorBehavior(ownedOrHalo))
    while isValid(iter)
        particle = Simulator.Iterators.:*(iter)
        if !insideLocalDomain(domain.localBoxMin, domain.localBoxMax, particle)
            push!(mp)
        end
        Simulator.Iterators.:++(iter)
    end

    # transfer particles
    pL = []
    pR = []
    if domain.rank == 0
        iter = AutoPasInterface.begin(autoPasContainer, IteratorBehavior(ownedOrHalo))
        while isValid(iter)
            particle = Simulator.Iterators.:*(iter)
            if getID(particle) == 10
                push!(pR, particle)
            end
            Simulator.Iterators.:++(iter)
        end
    else
        iter = AutoPasInterface.begin(autoPasContainer, IteratorBehavior(ownedOrHalo))
        while isValid(iter)
            particle = Simulator.Iterators.:*(iter)
            if getID(particle) == 11
                push!(pL, particle)
            end
            Simulator.Iterators.:++(iter)
        end
    end

    newParticles = sendAndReceiveParticles(pL, pR, domain.rank)
    
    for p in newParticles
        addParticle(autoPasContainer, p)
    end
end

function sendAndReceiveParticles(particlesLeft, particlesRight, rank)

    particlesLeftSerialized = serializeParticles(particlesLeft)
    particlesRightSerialized = serializeParticles(particlesRight)
    
    left, right = getLeftAndRightNeighbour(rank)

    leftReceiveSize, rightReceiveSize = exchangeBufferSizes(size(particlesLeftSerialized)[1], size(particlesRightSerialized)[1], left, right)
    
    leftReceive = Array{}(undef, leftReceive)
    rightReceive = Array{}(undef, rightReceive)

    rParticlesLeft = MPI.Irecv!(leftReceive, comm; source=left, tag=left+32) # TODO: check what tag is and maybe add it
    rParticlesRight = MPI.Irecv!(rightReceive, comm; source=right, tag=right+64) # TODO: check what tag is and maybe add it

    sParticleLeft = MPI.Isend(particlesLeftSerialized, comm; dest=left, tag=left+32) # TOOD: check what tag is and maybe add it 
    sParticleRight = MPI.Isend(particlesRightSerialized, comm; dest=right, tag=right+64) # TOOD: check what tag is and maybe add it 

    MPI.waitall([rParticlesLeft, rParticlesRight, sParticleLeft, sParticleRight])

    deserializedLeft = deserializeParticles(leftReceive)
    deserializedRight = deserializeParticles(rightReceive)

    return append!(deserializedLeft, deserializedRight)
end

function groupParticles(particles)
end

function determineMigratingParticles(autoPasContainer)
#
end

function exchangeBufferSizes(leftSize, rightSize, left, right)

    sendSizeLeft = Array{Int64}(undef, 1)
    sendSizeLeft[1] = leftSize
    sendSizeRight = Array{Int64}(undef, 1)
    sendSizeRight[1] = rightSize

    receiveSizeLeft = Array{Int64}(undef, 1)
    receiveSizeRight = Array{Int64}(undef, 1)

    rSizeLeft = MPI.Irecv!(receiveSizeLeft, comm; source=left, tag=left+32) # TODO: check what tag is and maybe add it
    rSizeRight = MPI.Irecv!(receiveSizeRight, comm; soruce=right, tag=right+64) # TOOD: check what tag is and maybe add it 

    sSizeLeft = MPI.Isend(sendSizeLeft, comm; dest=left, tag=left+32) # TOOD: check what tag is and maybe add it 
    sSizeRight = MPI.Isend(sendSizeRight, comm; dest=right, tag=right+64) # TOOD: check what tag is and maybe add it 

    MPI.waitall([rSizeLeft, rSizeRight, sSizeLeft, sSizeRight])

    return receiveSizeLeft[1], receiveSizeRight[1]
end

function getLeftAndRightNeighbour(rank)
    if rank == 0
        return 1,1

    else
        return 0,0
    end
end

function getHaloParticles()

end

function getOutOfDomainParticles()

end

