
module Simulator
# /mnt/c/Users/laura/Documents/BA_SH/testB/AutoPas/build/examples/md-julia/
module Particles
  using CxxWrap
  @wrapmodule(joinpath("../../../../build/examples/md-julia/","libjulia_bindings.so"), :define_module_particles)
  # @wrapmodule(joinpath("/dss/dsshome1/03/ge49sib2/bachelor_thesis/AutoPas/build/examples/md-julia/","libjulia_bindings.so"), :define_module_particles)

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
    setP,
    setV,
    setF,
    getPosition,
    getVelocity,
    getV,
    getF,
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
  # @wrapmodule(joinpath("/dss/dsshome1/03/ge49sib2/bachelor_thesis/AutoPas/build/examples/md-julia/","libjulia_bindings.so"), :define_module_options)

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
  lc_c08,
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
  # @wrapmodule(joinpath("/dss/dsshome1/03/ge49sib2/bachelor_thesis/AutoPas/build/examples/md-julia/","libjulia_bindings.so"), :define_module_iterators)

  function __init__()
    @initcxx
  end
export
  isValid
end

module Properties
  using CxxWrap
  @wrapmodule(joinpath("../../../../build/examples/md-julia/","libjulia_bindings.so"), :define_module_properties)
  # @wrapmodule(joinpath("/dss/dsshome1/03/ge49sib2/bachelor_thesis/AutoPas/build/examples/md-julia/","libjulia_bindings.so"), :define_module_properties)

  function __init__()
    @initcxx
  end

export
  getMass
end

module AutoPasInterface
# /dss/dsshome1/03/ge49sib2/bachelor_thesis/AutoPas2/AutoPas/examples/md-julia/src/julia-simulator/main2.jl
  using CxxWrap
  @wrapmodule(joinpath("../../../../build/examples/md-julia/","libjulia_bindings.so"), :define_module_autopas)
  # @wrapmodule(joinpath("/dss/dsshome1/03/ge49sib2/bachelor_thesis/AutoPas/build/examples/md-julia/","libjulia_bindings.so"), :define_module_autopas)

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
  generateCubeGrid,
  generateCubeGridJulia,
  generateUniformGridJulia

include("./ForceCalculation.jl")
export
  updateForces

include("./InputStructs.jl")
export
  InputParameters,
  CubeGridInput,
  Thermostat,
  UniformGrid

include("./Simulator.jl")
export
  startSimulation,
  printSimulation,
  simulate,
  startSimulationEx

include("./TimeDiscretization.jl")
export
  updatePositions,
  updateVelocities,
  updatePositionJM,
  convertAndSo,
  updatePo,
  updatePositionParallel,
  updateVelocitiesParallel

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

include("./FMol.jl")
export
  Molecule

end