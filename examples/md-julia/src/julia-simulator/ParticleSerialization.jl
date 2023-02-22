# include("./SimulatorModules.jl")

# using .Simulator.Particles

function serializeParticles(particles)
# create array of particles serialized
    particleArray = Array{Float64}(undef, 15 * length(particles))
    baseIndex = 1
    for particle in particles
        particleArray[baseIndex] = getID(particle)
        position = getPosition(particle)
        particleArray[baseIndex + 1] = position[1]
        particleArray[baseIndex + 2] = position[2]
        particleArray[baseIndex + 3] = position[3]
        velocity = getVelocity(particle)
        particleArray[baseIndex + 4] = velocity[1]
        particleArray[baseIndex + 5] = velocity[2]
        particleArray[baseIndex + 6] = velocity[3]
        force = getForce(particle)
        particleArray[baseIndex + 7] = force[1]
        particleArray[baseIndex + 8] = force[2]
        particleArray[baseIndex + 9] = force[3]
        oldForce = getOldF(particle)
        particleArray[baseIndex + 10] = oldForce[1]
        particleArray[baseIndex + 11] = oldForce[2]
        particleArray[baseIndex + 12] = oldForce[3]
        particleArray[baseIndex + 13] = getTypeId(particle)
        # particleArray[baseIndex + 14] = getOwnerShipState(particle)
        particleArray[baseIndex + 14] = 1 # TODO: change
        baseIndex += 15
    end

    return particleArray
end

function deserializeParticles(particleArray)

    particles = Vector{MoleculeJ{Float64}}()
    np = size(particleArray)[1]
    baseIndex = 1
    while baseIndex <= np
        
        particle = MoleculeJ{Float64}(
            [particleArray[baseIndex+1], particleArray[baseIndex+2], particleArray[baseIndex+3]],
            [particleArray[baseIndex+4], particleArray[baseIndex+5], particleArray[baseIndex+6]],
            convert(Int64, particleArray[baseIndex]),
            convert(Int64, particleArray[baseIndex+13]))
        
        setForce(particle, [particleArray[baseIndex+7], particleArray[baseIndex+8], particleArray[baseIndex+9]])
        setOldF(particle, [particleArray[baseIndex+10], particleArray[baseIndex+11], particleArray[baseIndex+12]])
        # setOwnershipState(particle, [particleArray[baseIndex + 14]])
        # TODO: set ownership state
        push!(particles, particle)
        baseIndex += 15
    end

    return particles
end
#=
p1 = MoleculeJ{Float64}([1.1, 1.2, 1.3], [2.1, 2.2, 2.3], 2, 3)
setForce(p1, [3.1, 3.2, 3.3])
setOldF(p1, [5.1, 5.2, 5.3])

p2 = MoleculeJ{Float64}([1.1, 1.2, 1.3], [2.1, 2.2, 2.3], 5, 6)
setForce(p2, [3.1, 3.2, 3.3])
setOldF(p2, [5.1, 5.2, 5.3])

result = serializeParticles([p1, p2])
println(result)

particles = deserializeParticles(result)
println(toString(particles[1]))
println(toString(particles[2]))
=#