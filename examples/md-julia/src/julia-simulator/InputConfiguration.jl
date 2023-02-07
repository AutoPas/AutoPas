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
    factorBrownianMotion::Float64
    CubeGridInput() = new()   
end

"""
struct containing parameters for Thermostat
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
    container::Vector{Options.ContainerOptionValue}
    verletRebuildFrequency ::Int64
    verletSkinRadiusPerTimestep::Float64
    verletClusterSize::Int64
    selectorStrategy::Options.SelectorStrategyOption
    dataLayout::Vector{Options.DataLayoutOptionValue}
    traversal::Vector{Options.TraversalOptionValue}
    tuningStrategy::Options.TuningStrategyOption
    mpiStrategy::Options.MPIStrategyOption
    tuningInterval::Int64
    tuningSamples::Int64
    tuningMaxEvidence::Int64
    functor::String  # functor option e.g Lennard-Jones (12-6) 7
    newton3::Vector{Options.Newton3OptionValue}
    cutoff::Float64
    boxMin::Vector{Float64}
    boxMax::Vector{Float64}
    cellSize::Vector{Float64}
    deltaT::Float64
    iterations::Int64
    globalForce::Vector{Float64}
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