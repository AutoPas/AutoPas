using StaticArrays
mutable struct Molecule
    pos::SVector{3,Float64}
    vel::SVector{3,Float64}
    f::SVector{3,Float64}
    of::SVector{3,Float64}
    mId::Int64
    tId::Int64
end