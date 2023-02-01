"""
structs setting input parameters
"""
module InputConfiguration
# include("./SimulatorModules.jl")
using ..Options
"""
abstract type for all shapes of particle clusters
"""
abstract type ParticleObjectInput end

"""
derived struct of ParticleObject defining a CubeGrid
"""
mutable struct CubeGridInput <: ParticleObjectInput
    particlesPerDimension::Vector{Int64}
    particleSpacing::Float64
    bottomLeftCorner::Vector{Float64}
    velocity::Vector{Float64}
    particleType::Int64
    particleEpsilon::Float64
    particleSigma::Float64
    particleMass::Float64
    CubeGridInput() = new()   
end

"""
struct contianing parameters for Thermostat
"""
mutable struct Thermostat
    initialTemperature::Float32
    targetTemperature::Float32
    deltaTemperature::Float32
    thermostatInterval::Int32
    addBrownianMotion::Bool
    Thermostat() = new()
end

"""
struct containing all parameters for the simulation
"""
mutable struct InputParameters
    container::Vector{String} # vector of contianer options -> parsing needed 1
    verletRebuildFrequency ::Int64
    verletSkinRadiusPerTimestep::Float64
    verletClusterSize::Int64
    selectorStrategy::String # 2
    dataLayout::Vector{String} # 3
    traversal::Vector{String} # 4
    tuningStrategy::String # 5
    mpiStrategy::String # 6
    tuningInterval::Int64
    tuningSamples::Int64
    tuningMaxEvidence::Int64
    functor::String  # functor option e.g Lennard-Jones (12-6) 7
    newton3::Vector{String} # 8
    cutoff::Float64
    boxMin::Vector{Float64}
    boxMax::Vector{Float64}
    cellSize::Vector{Float64}  # check which values can be used
    deltaT::Float64
    iterations::Int64
    periodicBoundaries::Bool
    objects::Vector{ParticleObjectInput}
    thermostat::Thermostat 
    logLevel::String  # log level maybe string # 9
    noFlops::Bool
    noEndConfig::Bool # what does this mean?
    noProgressBar::Bool # what does this mean
    vtkFilename::String
    vtkWriteFrequency::String
    InputParameters() = new()
end

function parseContainerOption(inputParameters)
    containerOptions = Vector{Options.ContainerOptionValue}([])
    for option in inputParameters.container
        if cmp(option, "directSum") == 0
            push!(containerOptions, Options.directSum)
        elseif cmp(option, "linkedCells") == 0
            push!(containerOptions, Options.linkedCells)
        elseif cmp(option, "linkedCellsReferences") == 0
            push!(containerOptions, Options.linkedCellsReferences)
        elseif cmp(option, "varVerletListsAsBuild") == 0
            push!(containerOptions, Options.varVerletListsAsBuild)
        elseif cmp(option, "verletClusterLists") == 0
            push!(containerOptions, Options.verletClusterLists)
        elseif cmp(option, "verletLists") == 0
            push!(containerOptions, Options.verletLists)
        elseif cmp(option, "verletListsCells") == 0
            push!(containerOptions, Options.verletListsCells)
        elseif cmp(option, "pairwiseVerletLists") == 0
            push!(containerOptions, Options.pairwiseVerletLists)
        elseif cmp(option, "octree") == 0
            push!(containerOptions, Options.octree)
        else
            throw(ArgumentError("argument needs to be of of enum type"))
        end        
    end
    return containerOptions
end


export
    CubeGridInput
    Thermostat
    InputParameters
end