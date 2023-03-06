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
                deleteParticle(autoPasContainer, particle)
            end
            Simulator.Iterators.:++(iter)
        end
    else
        iter = AutoPasInterface.begin(autoPasContainer, IteratorBehavior(ownedOrHalo))
        while isValid(iter)
            particle = Simulator.Iterators.:*(iter)
            if getID(particle) == 11
                push!(pL, particle)
                deleteParticle(autoPasContainer, particle)
            end
            Simulator.Iterators.:++(iter)
        end
    end

    newParticles = sendAndReceiveParticles(pL, pR, domain.rank)
    
    for p in newParticles
        addParticle(autoPasContainer, p)
    end
end

function sendAndReceiveParticles(particlesLeft, particlesRight, domain, comm)

    particlesLeftSerialized = serializeParticles(particlesLeft)
    particlesRightSerialized = serializeParticles(particlesRight)
    
    left, right = getLeftAndRightNeighbour(domain)

    leftReceiveSize, rightReceiveSize = exchangeBufferSizes(size(particlesLeftSerialized)[1], size(particlesRightSerialized)[1], left, right, domain.rank, comm)
    
    leftReceive = Array{Float64}(undef, leftReceiveSize)
    rightReceive = Array{Float64}(undef, rightReceiveSize)

    rParticlesLeft = MPI.Irecv!(leftReceive, comm; source=left, tag=left+32) # TODO: check what tag is and maybe add it
    rParticlesRight = MPI.Irecv!(rightReceive, comm; source=right, tag=right+32) # TODO: check what tag is and maybe add it

    sParticleLeft = MPI.Isend(particlesLeftSerialized, comm; dest=left, tag=domain.rank+32) # TOOD: check what tag is and maybe add it 
    sParticleRight = MPI.Isend(particlesRightSerialized, comm; dest=right, tag=domain.rank+32) # TOOD: check what tag is and maybe add it 

    MPI.Waitall([rParticlesLeft, rParticlesRight, sParticleLeft, sParticleRight])

    deserializedLeft = deserializeParticles(leftReceive)
    deserializedRight = deserializeParticles(rightReceive)

    return append!(deserializedLeft, deserializedRight)
end

function groupParticles(particles)
end

function determineMigratingParticles(autoPasContainer)
#
end

function exchangeBufferSizes(leftSize, rightSize, left, right, rank, comm)

    sendSizeLeft = Array{Int64}(undef, 1)
    sendSizeLeft[1] = leftSize
    sendSizeRight = Array{Int64}(undef, 1)
    sendSizeRight[1] = rightSize

    receiveSizeLeft = Array{Int64}(undef, 1)
    receiveSizeRight = Array{Int64}(undef, 1)

    rSizeLeft = MPI.Irecv!(receiveSizeLeft, comm; source=left, tag=left+32) # TODO: check what tag is and maybe add it
    rSizeRight = MPI.Irecv!(receiveSizeRight, comm; source=right, tag=right+32) # TOOD: check what tag is and maybe add it 

    sSizeLeft = MPI.Isend(sendSizeLeft, comm; dest=left, tag=rank+32) # TOOD: check what tag is and maybe add it 
    sSizeRight = MPI.Isend(sendSizeRight, comm; dest=right, tag=rank+32) # TOOD: check what tag is and maybe add it 

    MPI.Waitall([rSizeLeft, rSizeRight, sSizeLeft, sSizeRight])

    return receiveSizeLeft[1], receiveSizeRight[1]
end

function getLeftAndRightNeighbour(domain)
    left = mod(domain.rank-1,domain.size)
    right = mod(domain.rank+1, domain.size)

    return left, right
end

function getHaloParticles()

end

function getOutOfDomainParticles()

end

function getHaloAreaLeft(domain)
    # TODO: check if this calculation is correct; is it starting at localMin or localMin - cutoff?
    if domain.localBoxMin == domain.globalBoxMin
        # at globalMin boundary in x direction
        hlmin = copy(domain.globalBoxMin)
        hlmin[1] = domain.globalBoxMax[1]
        hlmin[2] -= domain.cutoff
        hlmin[3] -= domain.cutoff

        hlmax = copy(domain.globalBoxMax)
        hlmax = hlmax .+ domain.cutoff

        return hlmin, hlmax
    else
        # at local boundary in x direction
        hlmin = copy(domain.localBoxMin)
        hlmin[2] -= domain.cutoff
        hlmin[3] -= domain.cutoff

        hlmax = copy(domain.localBoxMax)
        hlmax[1] = domain.localBoxMin[1]
        hlmax = hlmax .+ domain.cutoff

        return hlmin, hlmax
    end
end

function getHaloAreaRight(domain)
    # TODO: check if this calculation is correct; is it starting at localMin or localMin - cutoff?
    if domain.localBoxMax == domain.globalBoxMax
        # at globalMax boundary in x direction
        hrmin = domain.globalBoxMin .- domain.cutoff
        
        hrmax = copy(domain.globalBoxMax)
        hrmax[1] = domain.globalBoxMin[1]
        hrmax[2] += domain.cutoff
        hrmax[3] += domain.cutoff

        return hrmin, hrmax
    else
        # at local boundary in x direction
        hrmin = copy(domain.localBoxMin)
        hrmin[1] = domain.localBoxMax[1]
        hrmin = hrmin .- domain.cutoff

        hrmax = copy(domain.localBoxMax)
        hrmax[2] += domain.cutoff
        hrmax[3] += domain.cutoff

        return hrmin, hrmax
   
    end
end

"""
send and receive halo particles
"""

function exchangeHaloParticles(autoPasContainer, domain, comm)
    println("in exchange halo particles")
    # 0. check if at border of global domain and see what happens in this case ?
    # println("halo")
    # 1. iterate over region and collect the particles
    particlesLeft = []
    lMin, lMax = getHaloAreaLeft(domain)
    iter = regionIterator(autoPasContainer, lMin, lMax, IteratorBehavior(ownedOrHalo))
    while isValid(iter)
        println("in left halo")
        push!(particlesLeft, Simulator.Iterators.:*(iter))
        Simulator.Iterators.:++(iter)
    end

    particlesRight = []
    rMin, rMax = getHaloAreaRight(domain)
    iter = regionIterator(autoPasContainer, rMin, rMax, IteratorBehavior(ownedOrHalo))
    while isValid(iter)
        println("in right halo")
        push!(particlesRight, Simulator.Iterators.:*(iter))
        Simulator.Iterators.:++(iter)
    end

    # 2. send and receive halo particles

    incomingHalo = sendAndReceiveParticles(particlesLeft, particlesRight, domain, comm)

    # 3. add particles to container (with addHalo -> check if addHalo is different to other methods)
    for particle in incomingHalo

        # TODO: check if addhalo does the same
        # addHaloParticle(autoPasContainer, particle)
        setOwnershipState(particle, haloS)
        addHaloParticle(autoPasContainer, particle)
    end

end

"""
send and receive particles migrated from rank domain i to rank domain j with i != j
"""

function exchangeMigratingParticles(autoPasContainer, domain, comm)
    println("in exchange migrating particles")
    # calculate migrating area for left neighbour
    mLeftMin = copy(domain.localBoxMin)
    mLeftMin = mLeftMin .- domain.cutoff
    mLeftMax = copy(domain.localBoxMax)
    mLeftMax[1] -= domain.rankWidth
    mLeftMax[2] += domain.cutoff
    mLeftMax[3] += domain.cutoff

    
    # 1. iterate over region and collect the particles
    particlesLeft = []
    iter = regionIterator(autoPasContainer, mLeftMin, mLeftMax, IteratorBehavior(ownedOrHalo))
    while isValid(iter)
        println("in left area")
        particle = Simulator.Iterators.:*(iter)
        pos = getPosition(particle)
        if mLeftMin < domain.globalBoxMin
            for i in 1:3
                if pos[i] < domain.globalBoxMin[i]
                    pos[i] += domain.dimensions[i]
                end
            end
            setPosition(particle, pos)
        end
        push!(particlesLeft, particle)
        # deleteParticle(autoPasContainer, particle)
        Simulator.Iterators.:++(iter)
    end

    # calculate migrating area for right neighbour
    #=
    mRightMin = copy(domain.localBoxMin)
    mRightMin[1] += domain.rankWidth
    mRightMin[2] -= domain.cutoff
    mRightMin[3] -= domain.cutoff
    mRightMax = copy(domain.localBoxMax)
    mRightMax = mRightMax .+ domain.cutoff
    =#
    mRightMin = mLeftMin
    mRightMin[1] += (domain.rankWidth + domain.cutoff)
    mRightMax = mLeftMax
    mRightMax[1] += (domain.rankWidth + domain.cutoff)

    particlesRight = []
    iter = regionIterator(autoPasContainer, mRightMin, mRightMax, IteratorBehavior(ownedOrHalo))
    while isValid(iter)
        println("in right area")
        # TODO: if at boundary -> change position of particle
        particle = Simulator.Iterators.:*(iter)
        pos = getPosition(particle)
        
        if pos > domain.globalBoxMax
            for i in 1:3
                if pos[i] > domain.globalBoxMax[i]
                    pos[i] = pos[i] - domain.dimensions[i]
                end
            end
            setPosition(particle, pos)
        end
        push!(particlesRight, particle)
        # deleteParticle(autoPasContainer, particle)
        Simulator.Iterators.:++(iter)
    end

    # 2. send and receive migrating particles

    incoming = sendAndReceiveParticles(particlesLeft, particlesRight, domain, comm)
    
    for p in particlesLeft
        deleteParticle(autoPasContainer, p)
    end

    for p in particlesRight
        deleteParticle(autoPasContainer, p)
    end

    # 3. add particles to container
    for particle in incoming
        addParticle(autoPasContainer, particle)
    end
    
end

function moveParticles(autoPasContainer, domain, globalBorder, localBorder)
    cutoff = 1.5
    for dimension in 1:3
        if (( dimension == 1) && (gloablBorder == localBorder) ) || (dimension != 1)
            
            shiftArray = zeros(3)
            shiftArray[dimension] = domain.dimensions[dimension]

            lowerMinCorner = zeros(3)
            lowerMaxCorner = zeros(3)
            upperMinCorner = zeros(3)
            upperMaxCorner = zeros(3)

            for i in 1:3
                upperMaxCorner[i] = domain.globalBoxMax[i] + cutoff
                lowerMinCorner[i] = domain.globalBoxMin[i] - cutoff
                if i != dimension
                    upperMinCorner[i] = domain.globalBoxMin[i] - cutoff
                    lowerMaxCorner[i] = domain.globalBoxMax[i] + cutoff
                else
                    upperMinCorner[i] = domain.globalBoxMax[i]
                    lowerMaxCorner[i] = domain.globalBoxMin[i]
                end
            end
        end
    end
end

function moveOuterParticles(autoPasContainer, domain)
    println("in move outer particles")
    # create outofbox domains for y and z direction
    cutoff = 1.5
    li = domain.localBoxMin
    la = domain.localBoxMax
    for dimension in 1:3
        shiftArray = [0.0, 0.0, 0.0]
        shiftArray[dimension] = domain.dimensions[dimension]
        upperMinCorner = [0.0, 0.0, 0.0]
        upperMaxCorner = [0.0, 0.0, 0.0]

        lowerMinCorner = [0.0, 0.0, 0.0]
        lowerMaxCorner = [0.0, 0.0, 0.0]

        for j in 1:3
            upperMaxCorner[j] = la[j] + domain.cutoff
            lowerMinCorner[j] = li[j] - domain.cutoff
            if j != dimension
                upperMinCorner[j] = li[j] - domain.cutoff
                lowerMaxCorner[j] = la[j] + domain.cutoff
            else
                upperMinCorner[j] = la[j]
                lowerMaxCorner[j] = li[j]
            end
        end

        if ((dimension == 1) && (domain.localBoxMin == domain.globalBoxMin)) || (dimension != 1)  
            iter = regionIterator(autoPasContainer, lowerMinCorner, lowerMaxCorner, IteratorBehavior(ownedOrHalo))

            while isValid(iter)
                println("###########moved particle lower\n")
                # particle = Simulator.Iterators.:*(iter)
                addPosition(Simulator.Iterators.:*(iter), shiftArray)
                Simulator.Iterators.:++(iter)
            end 
        end
        if ((dimension == 1) && (domain.localBoxMax == domain.globalBoxMax)) || (dimension != 1)  
            iter = regionIterator(autoPasContainer, upperMinCorner, upperMaxCorner, IteratorBehavior(ownedOrHalo))

            while isValid(iter)
                println("###############moved particle upper\n")
                # particle = Simulator.Iterators.:*(iter)
                addPosition(Simulator.Iterators.:*(iter), (-1) * shiftArray)
                Simulator.Iterators.:++(iter)
            end
        end   
    end
    return autoPasContainer

end
 
function handleXBoundary(autoPasContainer, domain, particle, index)

    tmp = zeros(3)
    pos = getPosition(particle)
    isCorner = 0
    if pos[1] < (domain.globalBoxMin[1] + domain.cutoff)
        tmp[1] += domain.dimensions[1]
        isCorner += 1

    elseif pos[1] > (domain.globalBoxMax[1] - domain.cutoff)
        tmp[1] -= domain.dimensions[1]
        isCorner += 1
    end
    
    if pos[2] < (domain.globalBoxMin[2] + domain.cutoff)
        tmp[2] += domain.dimensions[2]
        isCorner += 1
        xyPos = copy(pos)
        xyPos = xyPos + tmp
        xyHalo = MoleculeJ{Float64}(xyPos, getVelocity(particle), index, getTypeId(particle))
        index += 1
        addHaloParticle(autoPasContainer, xyHalo)

    elseif pos[2] > (domain.globalBoxMax[2] - domain.cutoff)
        tmp[2] -= domain.dimensions[2]
        isCorner += 1
        xyPos = copy(pos)
        xyPos = xyPos + tmp
        xyHalo = MoleculeJ{Float64}(xyPos, getVelocity(particle), index, getTypeId(particle))
        index += 1
        addHaloParticle(autoPasContainer, xyHalo)
    end

    if pos[3] < (domain.globalBoxMin[3] + domain.cutoff)
        tmp[3] += domain.dimensions[3]
        isCorner += 1
        xzPos = copy(pos)
        xzPos[1] += tmp[1]
        xzPos[3] += tmp[3]
        xzHalo = MoleculeJ{Float64}(xzPos, getVelocity(particle), index, getTypeId(particle))
        index += 1
        addHaloParticle(autoPasContainer, xzHalo)

    elseif pos[3] > (domain.globalBoxMax[3] - domain.cutoff)
        tmp[3] -= domain.dimensions[3]
        isCorner += 1
        xzPos = copy(pos)
        xzPos[1] += tmp[1]
        xzPos[3] += tmp[3]
        xzHalo = MoleculeJ{Float64}(xzPos, getVelocity(particle), index, getTypeId(particle))
        index += 1
        addHaloParticle(autoPasContainer, xzHalo)
    end
    
    if isCorner == 3
        xyzPos = copy(pos)
        xyzPos = xyzPos + tmp
        xyzHalo = MoleculeJ{Float64}(xyzPos, getVelocity(particle), index, getTypeId(particle))
        index += 1
        addHaloParticle(autoPasContainer, xyzHalo)
    end
    return index
end

"""
insert halo particles from particles in upper boundary zone
"""

function insertHaloParticles(autoPasContainer, minCorner, maxCorner, i, domain, index)
    println("in insertHaloParticlesUpper")
    # iter = regionIterator(autoPasContainer, minCorner, maxCorner, IteratorBehavior(owned))
    iter = Simulator.AutoPasInterface.begin(autoPasContainer, IteratorBehavior(ownedOrHalo))
    np = getNumberOfParticles(autoPasContainer, IteratorBehavior(ownedOrHalo))
    println("number of partiles: $np")
    # println("$minCorner; $maxCorner")
    wasInLoop = false
    while isValid(iter)
        wasInLoop = true
        println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        # println(toString(Simulator.Iterators.:*(iter)))
        # add halo particle at dimension i
        particle = Simulator.Iterators.:*(iter)
        ps = toString(particle)
        rank = domain.rank
        println("$rank: $ps")
        pos = getPosition(particle)
        newP = copy(pos)
        if newP[i] > (domain.globalBoxMax[i] - domain.cutoff)
            newP[i] -= domain.dimensions[i]
        else
            newP[i] += domain.dimensions[i]
        end

        println("new pos: $newP")
       
        pHalo = MoleculeJ{Float64}(newP , getVelocity(particle), index, getTypeId(particle))
        index += 1
        # addHaloParticle(autoPasContainer, pHalo)
        println("added halo")
        
        if i == 1
            index = handleXBoundary(autoPasContainer, domain, particle, index)
        end
        
        if i == 2
            # check if particle is in edge of y-z plane
            yzPos = copy(newP)
            isOutOfZ = false
            if pos[3] < (domain.globalBoxMin[3] + domain.cutoff)
                isOutOfZ = true
                yzPos[3] += domain.dimensions[3]
            elseif pos[3] > (domain.globalBoxMax[3] - domain.cutoff)
                isOutOfZ = true
                yzPos[3] += domain.dimensions[3]
            end
            
            if isOutOfZ
                yzHalo = MoleculeJ{Float64}(yzPos, getVelocity(particle), index, getTypeId(particle))
                index += 1
                addHaloParticle(autoPasContainer, yzHalo)
            end

        end
        Simulator.Iterators.:++(iter)
    end
    println("was in loop: ", wasInLoop)
    return index
end

"""
insert halo particles from particles in lower boundary zone
"""

function insertLowerHaloParticles(autoPasContainer, minCorner, maxCorner, i, domain, index)
    println("in insertHaloParticlesLower")
    # iter = Simulator.AutoPasInterface.begin(autoPasContainer, IteratorBehavior(ownedOrHalo)) # regionIterator(autoPasContainer, minCorner, maxCorner, IteratorBehavior(ownedOrHalo))
    iter = regionIterator(autoPasContainer, minCorner, maxCorner, IteratorBehavior(ownedOrHalo))
    # println("$minCorner; $maxCorner")
    wasInLoop = false
    while isValid(iter)
        wasInLoop = true
        println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        # println(toString(Simulator.Iterators.:*(iter)))
        # add halo particle at dimension i
        particle = Simulator.Iterators.:*(iter)
        pos = getPosition(particle)
        newP = copy(pos)
        newP[i] += domain.dimensions[i] # mod(newP[i] + domain.dimensions[i], domain.dimensions[i])

        println("new pos: $newP")
       
        pHalo = MoleculeJ{Float64}(newP , getVelocity(particle), convert(UInt, index), getTypeId(particle))
        index += 1
        addHaloParticle(autoPasContainer, pHalo)
        println("added halo")
        
        if i == 1
            handleXBoundary(autoPasContainer, domain, particle)
        end
        
        #=
        if i == 2
            # check if particle is in edge of y-z plane
            if (pos[2] < (li[2] + cutoff)) && (pos[3] < (li[3] + cutoff))
                newPos = zeros(3)
                tmp = [0.0, domain.dimensions[2], domain.dimensions[3]]
                for i in 1:3
                    newPos[i] = mod(pos[i] + tmp[i], domain.dimensions[i])
                end
                pHaloExz = MoleculeJ{Float64}(newPos, getVelocity(particle), getID(particle), getTypeId(particle))
                addHaloParticle(autoPasContainer, pHaloExz)
            end

        end
        =#
        Simulator.Iterators.:++(iter)
    end
    println("was in loop: ", wasInLoop)
    return index
end

function insertHaloParticles(autoPasContainer, domain, index)
    println("in insertHaloParticles")
    cutoff = 1.5
    li = domain.localBoxMin
    la = domain.localBoxMax
    
    for i in 1:3
        if domain.localBoxMin == domain.globalBoxMin
        end
        shiftArray = [0.0, 0.0, 0.0]
        shiftArray[i] = domain.dimensions[i]

        # set domain sizes
        lowerMinCorner = copy(li)
        lowerMaxCorner = copy(la)
        lowerMaxCorner[i] = li[i] + domain.cutoff

        println("lower: $lowerMinCorner ; $lowerMaxCorner")

        upperMinCorner = copy(li)
        upperMinCorner[i] = la[i] - domain.cutoff
        upperMaxCorner = copy(la)

        println("upper: $upperMinCorner ; $upperMaxCorner")
        if ((i == 1) && (domain.localBoxMax == domain.globalBoxMax)) || (i != 1)
            println("in upper")
            index = insertHaloParticles(autoPasContainer, upperMinCorner, upperMaxCorner, i, domain, index)
        end

        if ((i == 1) && (domain.localBoxMin == domain.globalBoxMin)) || (i != 1)
            println("in lower")
            index = insertHaloParticles(autoPasContainer, lowerMinCorner, lowerMaxCorner, i, domain, index)
        end
        
    end
    return index
end


"""
apply periodic boundary condition on all rank boundaries
"""

function applyPeriodicBoundary(autoPasContainer, domain, index)
    println("in applyBoundary condition")
    autoPasContainer = moveOuterParticles(autoPasContainer, domain)

    println("starting particle print: ")
    iter = Simulator.AutoPasInterface.begin(autoPasContainer, Options.IteratorBehavior(Options.ownedOrHalo))

    while isValid(iter)
        particle = Simulator.Iterators.:*(iter)
        # ost = ost*toString(particle)*"\n"
        println(toString(particle))
        Simulator.Iterators.:++(iter)
    end
    println("ending particle print\n")

    in = insertHaloParticles(autoPasContainer, domain, index)
    return in
end

function deleteHaloParticlesLocalBounds(autoPasContainer, domain)
    iter = regionIterator(autoPasContainer, domain.mLeftMin, domain.mLeftMax, IteratorBehavior(ownedOrHalo))

    while isValid(iter)
        particle = Simulator.Iterators.:*(iter)
        if getOwnershipState(particle) == haloS
            deleteParticle(autoPasContainer, particle)
        end
        Simulator.Iterators.:++(iter)
    end

    iter = regionIterator(autoPasContainer, domain.mRightMin, domain.mRightMax, IteratorBehavior(ownedOrHalo))

    while isValid(iter)
        particle = Simulator.Iterators.:*(iter)
        if getOwnershipState(particle) == haloS
            deleteParticle(autoPasContainer, particle)
        end
        Simulator.Iterators.:++(iter)
    end
end

