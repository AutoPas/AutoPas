include("./Domain.jl")
include("./DomainDecomposition.jl")
using .Domain, .DomainDecomposition


function getHaloAreaLeft(domain)
    # TODO: check if this calculation is correct; is it starting at localMin or localMin - cutoff?
    if domain.localBoxMin == domain.globalBoxMin
        # at globalMin boundary in x direction
        hlmin = copy(domain.globalBoxMin)
        hlmin[1] = globalBoxMax[1]
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
        hrmin = globalBoxMin .- domain.cutoff
        
        hrmax = copy(globalBoxMax)
        hrmax[1] = globalBoxMin[1]
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

function getOutOfBounds(domain)
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
        println("$dimension: $lowerMinCorner ; $lowerMaxCorner")
        println("$dimension: $upperMinCorner ; $upperMaxCorner")
    end
end
globalBoxMin = [0.0, 0.0, 0.0]
globalBoxMax = [10.0, 10.0, 10.0]


localBoxMin = [0.0, 0.0, 0.0]
localBoxMax = [5.0, 10.0, 10.0]

width = localBoxMax[1] - localBoxMin[1]

dimensions = globalBoxMax - globalBoxMin

domain = Domain(0, 2, width, 1.5, localBoxMin, localBoxMax, globalBoxMin, globalBoxMax, dimensions)

#= 
lMin, lMax = getHaloAreaLeft(domain)
println("lmin: $lMin")
println("lmax: $lMax")

rMin, rMax = getHaloAreaRight(domain)

println("rmin: $rMin")
println("rmax: $rMax")
=#

getOutOfBounds(domain)