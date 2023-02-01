module Generator

using ..Particles
using ..InputConfiguration


function generateCubeGrid(cubeGrid)
    particles = Vector{Particles.MoleculeJ{Float64}}([])
    id_ = 0
    for x in 1 : cubeGrid.particlesPerDimension[1]
        for y in 1 : cubeGrid.particlesPerDimension[2]
            for z in 1 : cubeGrid.particlesPerDimension[3]
                
                pos = [(x-1) * cubeGrid.particleSpacing, (y-1) * cubeGrid.particleSpacing, (z-1) * cubeGrid.particleSpacing] + cubeGrid.bottomLeftCorner
                
                particle = Particles.MoleculeJ{Float64}()
                Particles.setPos(particle, pos)
                push!(particles, particle)  

            end
        end
    end
    return particles
end

export
    generateCubeGrid
end