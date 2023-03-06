using MPI
MPI.Init()

mutable struct Domain
    rank::Int64
    size::Int64
    rankWidth::Float64
    cutoff::Float64
    localBoxMin::Vector{Float64}
    localBoxMax::Vector{Float64}
    globalBoxMin::Vector{Float64}
    globalBoxMax::Vector{Float64}
    #=
    hLeftMin::Vector{Float64} # "halo" cell of left neighbor min corner
    hLeftMax::Vector{Float64} # "halo" cell of left neighbor max corner
    hRightMin::Vector{Float64} # "halo" cell of right neighbor min corner
    hRightMax::Vector{Float64} # "halo" cell of right neighbor max corner
    mLeftMin::Vector{Float64} # migrating cell of left neighbor min corner
    mLeftMax::Vector{Float64} # migrating cell of left neighbor max corner
    mRightMin::Vector{Float64} # migrating cell of right neighbor min corner
    mRightMax::Vector{Float64} # migrating cell of right neighbor max corner
    =#
    dimensions::Vector{Float64} # length of x, y and z dimension, where x == rankwidth

end

function createDomain(globalBoxMin, globalBoxMax, rankWidth, comm, cutoff)
    rankId = MPI.Comm_rank(comm)# MPI.Comm_rank(comm)
    
    lMin = globalBoxMin[1] + rankId * rankWidth
    lMax = globalBoxMin[1] + (rankId+1) * rankWidth
    localBoxMin = [globalBoxMin[1] + rankId * rankWidth, globalBoxMin[2], globalBoxMin[3]]
    localBoxMax = [globalBoxMin[1] + (rankId+1) * rankWidth, globalBoxMax[2], globalBoxMax[3]]

    dx = localBoxMax[1] - localBoxMin[1]

    tmp = [0.0, 0.0, 0.0]
    
    hlmin = [localBoxMin[1], localBoxMin[2], localBoxMin[3]]
    hlmax = [localBoxMin[1] + cutoff, localBoxMax[2], localBoxMax[3]]

    hrmin = [localBoxMax[1] - cutoff, localBoxMin[2], localBoxMin[3]]
    hrmax = [localBoxMax[1], localBoxMax[2], localBoxMax[3]]

    # println("boxMin: $globalBoxMin")
    # println("box Max: $globalBoxMax")
    mlmin = [localBoxMin[1] - cutoff, localBoxMin[2], localBoxMin[3]]
    mlmax = [localBoxMax[1] - dx, localBoxMax[2], localBoxMax[3]]

    mrmin = [localBoxMin[1] + dx, localBoxMin[2], localBoxMin[3]]
    mrmax = [localBoxMax[1] + cutoff, localBoxMax[2], localBoxMax[3]]

    # TODO: add globalBoxMin/Max

    dimensions = globalBoxMax - globalBoxMin

    domain = Domain(rankId, MPI.Comm_size(comm), rankWidth, cutoff, localBoxMin, localBoxMax, globalBoxMin, globalBoxMax, dimensions)
    # Domain(rankId, MPI.Comm_size(comm), rankWidth, localBoxMin, localBoxMax, globalBoxMin, globalBoxMax, hlmin, hlmax, hrmin, hrmax, mlmin, mlmax, mrmin, mrmax, tmp)

    return domain
end

function calculateMPIWidth(globalBoxMin, globalBoxMax, comm)
    return (globalBoxMax[1] - globalBoxMin[1]) / MPI.Comm_size(comm) 
end

function insideLocalDomain(localBoxMin, localBoxMax, particle)
    pos = getPosition(particle)
    for i in 1:3
        if (pos[i] < localBoxMin[i]) || (pos[i] > localBoxMax[1])
            return false
        end
    end
    return true
end

#=
b1 = [0.0, 0.0, 0.0]
b2 = [10.0, 10.0, 10.0]
comm = MPI.COMM_WORLD

width = calculateMPIWidth(b1, b2, comm)

sb = createDomain(b1, b2, width, comm, 1.0)
=#
# println("rank: $sb \n")
# println(sb.localBoxMin)
# println(sb.localBoxMax)