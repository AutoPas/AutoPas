using Random, Distributions

function generateCubeGrid(cubeGrid, mId)
    particles = Vector{Particles.MoleculeJ{Float64}}([])
    id_ = 0
    generator = createRandomGenerator()
    println("generation")
    for x in 1 : cubeGrid.particlesPerDimension[1]
        for y in 1 : cubeGrid.particlesPerDimension[2]
            for z in 1 : cubeGrid.particlesPerDimension[3]
                pos = [(x-1) * cubeGrid.particleSpacing, (y-1) * cubeGrid.particleSpacing, (z-1) * cubeGrid.particleSpacing] + cubeGrid.bottomLeftCorner
                println(pos)
                # particle = Particles.MoleculeJ{Float64}()
                # Particles.setPos(particle, pos)
                v = cubeGrid.velocity + addBrownianMotion(cubeGrid.factorBrownianMotion, generator)
                println(v)
                particle = Particles.MoleculeJ{Float64}(pos, v, mId, cubeGrid.particleType)
                # particle = Particles.MoleculeJ{Float64}()
                # Particles.setPos(particle, pos)
                # Particles.setV(particle, v)
                println(Particles.getPos(particle))
                println(Particles.getV(particle))
                # println(Particles.toString(particle))
                push!(particles, particle)
                mId = mId + 1
                println("########")
            end
        end
    end
    println("done generation")
    for p in particles
        # println(Particles.toString(p))
    end
    return particles, mId
end

function generateObject(particleObject, mId)

    if typeof(particleObject) == typeof(CubeGridInput())
        return generateCubeGrid(particleObject, mId)
    end

end

function createRandomGenerator()

    Random.seed!(35)
    return Normal()

end

function addBrownianMotion(average_velocity, randomGenerator)::Vector{Float64}
    return average_velocity * rand(randomGenerator, 3)
end