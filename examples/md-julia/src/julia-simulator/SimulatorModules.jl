
module Simulator

module Particles
  using CxxWrap
  @wrapmodule(joinpath("../../../../build/examples/md-julia/","libjulia_bindings.so"), :define_module_particles)

  function __init__()
    @initcxx
  end

  export
    MoleculeJ,
    setPosition,
    setVelocity,
    setForce,
    setOldF,
    setID,
    setTypeId,
    getPosition,
    getVelocity,
    getForce,
    getOldF,
    getID,
    getTypeId,
    addPosition,
    addVelocity,
    addForce,
    subForce,
    toString,
    setOwnershipState,
    getOwnershipState,
    ownedS,
    haloS

end

module Options
  using CxxWrap
  @wrapmodule(joinpath("../../../../build/examples/md-julia/","libjulia_bindings.so"), :define_module_options)

  function __init__()
    @initcxx
  end
export
  ContainerOption,
  DataLayoutOption,
  Newton3Option,
  TraversalOption,
  SelectorStrategyOption,
  TuningStrategyOption,
  MPIStrategyOption,
  IteratorBehavior
export
  linkedCells,
  directSum,
  fastestAbs,
  aos,
  lc_c01,
  fullSearch,
  noMPI,
  disabled,
  enabled,
  owned,
  halo,
  ownedOrHalo
end

module Iterators
  using CxxWrap
  @wrapmodule(joinpath("../../../../build/examples/md-julia/","libjulia_bindings.so"), :define_module_iterators)

  function __init__()
    @initcxx
  end
export
  isValid
end

module Properties
  using CxxWrap
  @wrapmodule(joinpath("../../../../build/examples/md-julia/","libjulia_bindings.so"), :define_module_properties)

  function __init__()
    @initcxx
  end

export
  getMass
end

module AutoPasInterface
  using CxxWrap
  @wrapmodule(joinpath("../../../../build/examples/md-julia/","libjulia_bindings.so"), :define_module_autopas)

  function __init__()
    @initcxx
  end
  export
    updateContainer,
    AutoPas,
    init,
    setAllowedContainers,
    setAllowedDataLayouts,
    setAllowedNewton3Options,
    setAllowedTraversals,
    setBoxMin,
    setBoxMax,
    setCutoff,
    setSelectorStrategy,
    setTuningInterval,
    setTuningStrategyOption,
    setMPIStrategy,
    addParticle,
    iteratePairwise,
    getNumberOfParticles,
    updateContainer,
    getCutoff,
    deleteParticle,
    getBoxMin,
    getBoxMax,
    regionIterator,
    addHaloParticle
end

using .Particles, .Options, .Iterators, .Properties, .AutoPasInterface
include("./Configuration.jl")
export
  parseInput,
  getNp

include("./Generator.jl")
export
  generateObject,
  generateCubeGrid

include("./ForceCalculation.jl")
export
  updateForces

include("./InputConfiguration.jl")
export
  InputParameters,
  CubeGridInput,
  Thermostat

include("./Simulator.jl")
export
  startSimulation,
  printSimulation,
  simulate

include("./TimeDiscretization.jl")
export
  updatePositions,
  updateVelocities

include("./Domain.jl")
export
  Domain,
  createDomain,
  calculateMPIWidth,
  insideLocalDomain

include("./DomainDecomposition.jl")
export
  migrateParticles,
  sendAndReceiveParticles,
  groupParticles,
  determineMigratingParticles,
  exchangeBufferSizes,
  getLeftAndRightNeighbour,
  exchangeHaloParticles,
  exchangeMigratingParticles,
  deleteHaloParticlesLocalBounds,
  applyPeriodicBoundary

include("./ParticleSerialization.jl")
export
  serializeParticles,
  deserializeParticles

end