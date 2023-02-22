using MPI
MPI.Init()

mutable struct Domain
    rank::Int64
    rankWidht::Float64
    localBoxMin::Vector{Float64}
    localBoxMax::Vector{Float64}
end

function createDomain(globalBoxMin, globalBoxMax, rankWidth, comm)
    rankId = MPI.Comm_rank(comm)

    lMin = globalBoxMin[1] + rankId * rankWidth
    lMax = globalBoxMin[1] + (rankId+1) * rankWidth
    
    globalBoxMin[1] = lMin
    globalBoxMax[1] = lMax 

    domain = Domain(rankId, rankWidth, globalBoxMin, globalBoxMax)

    return domain
end

function calculateMPIWidth(globalBoxMin, globalBoxMax, comm)
    return (globalBoxMax[1] - globalBoxMin[1]) / MPI.Comm_size(comm) 
end

function insideLocalDomain(localBoxMin, localBoxMax, particle)
    return false
end

b1 = [0.0, 0.0, 0.0]
b2 = [10.0, 10.0, 10.0]
comm = MPI.COMM_WORLD

width = calculateMPIWidth(b1, b2, comm)

sb = createDomain(b1, b2, width, comm)

println("rank: $sb \n")
# println(sb.localBoxMin)
# println(sb.localBoxMax)