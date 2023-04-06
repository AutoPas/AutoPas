using Random, Distributions
# using ..Simulator.Particles

function generateCubeGrid(cubeGrid, mId)
    particles = Vector{MoleculeJ{Float64}}([])
    id_ = 0
    generator = createRandomGenerator()
    for z in 1 : cubeGrid.particlesPerDimension[3]
        for y in 1 : cubeGrid.particlesPerDimension[2]
            for x in 1 : cubeGrid.particlesPerDimension[1]
                pos = [(x-1) * cubeGrid.particleSpacing, (y-1) * cubeGrid.particleSpacing, (z-1) * cubeGrid.particleSpacing] + cubeGrid.bottomLeftCorner
                v = cubeGrid.velocity # + addBrownianMotion(cubeGrid.factorBrownianMotion, generator)
                particle = MoleculeJ{Float64}(pos, v, mId, cubeGrid.particleType)
                push!(particles, particle)
                mId = mId + 1
            end
        end
    end
    return particles, mId
end

function generateCubeGridJulia(cubeGrid, mId)
    if typeof(cubeGrid) != typeof(CubeGridInput())
        return [], mId
    end

    particles = Vector{Molecule}([])
    id_ = 0
    generator = createRandomGenerator()
    for x in 1 : cubeGrid.particlesPerDimension[1]
        for y in 1 : cubeGrid.particlesPerDimension[2]
            for z in 1 : cubeGrid.particlesPerDimension[3]
                pos = [(x-1) * cubeGrid.particleSpacing, (y-1) * cubeGrid.particleSpacing, (z-1) * cubeGrid.particleSpacing] + cubeGrid.bottomLeftCorner
                v = cubeGrid.velocity + addBrownianMotion(cubeGrid.factorBrownianMotion, generator)
                particle = Molecule(pos, v, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], mId, cubeGrid.particleType)
                push!(particles, particle)
                mId = mId + 1
            end
        end
    end
    return particles, mId
end

function generateUniformGridJulia(uniformGrid, mId)
    particles = Vector{Molecule}([])
    id_ = 0
    generator = createRandomGenerator()
    for i in 1:uniformGrid.numParticles
        distVec = getUniformDistribution()
        pos = uniformGrid.bottomLeftCorner + uniformGrid.boxLength .* distVec
        v = uniformGrid.velocity + addBrownianMotion(uniformGrid.factorBrownianMotion, generator)
        particle = Molecule(pos, v, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], mId, uniformGrid.particleType)
        push!(particles, particle)
        mId += 1
    end
    return particles, mId
end

function generateUniformGrid(uniformGrid, mId)
    particles = Vector{MoleculeJ{Float64}}([])
    id_ = 0
    generator = createRandomGenerator()
    for i in 1:uniformGrid.numParticles
        distVec = getUniformDistribution()
        pos = uniformGrid.bottomLeftCorner + uniformGrid.boxLength .* distVec
        v = uniformGrid.velocity + addBrownianMotion(uniformGrid.factorBrownianMotion, generator)
        particle = MoleculeJ{Float64}(pos, v, mId, uniformGrid.particleType)
        push!(particles, particle)
        mId += 1
    end
    return particles, mId
end

function generateObject(particleObject, mId)

    if typeof(particleObject) == typeof(CubeGridInput())
        return generateCubeGrid(particleObject, mId)
    end

    if typeof(particleObject) == typeof(UniformGrid())
        return generateUniformGrid(particleObject, mId)
    end

end

function createRandomGenerator()
    Random.seed!(35)
    return Normal()

end

function addBrownianMotion(average_velocity, randomGenerator)::Vector{Float64}
    return average_velocity * rand(randomGenerator, 3)
end

function getUniformDistribution()
    return rand(3)
end