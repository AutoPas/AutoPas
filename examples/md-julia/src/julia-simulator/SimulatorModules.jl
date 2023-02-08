
module Particles
  using CxxWrap
  @wrapmodule(joinpath("../../../../build/examples/md-julia/","libjulia_bindings.so"), :define_module_particles)

  function __init__()
    @initcxx
  end
  export
    MoleculeJ
end

module Options
  using CxxWrap
  @wrapmodule(joinpath("../../../../build/examples/md-julia/","libjulia_bindings.so"), :define_module_options)

  function __init__()
    @initcxx
  end
end

module AIterators
  using CxxWrap
  @wrapmodule(joinpath("../../../../build/examples/md-julia/","libjulia_bindings.so"), :define_module_iterators)

  function __init__()
    @initcxx
  end
end

module Properties
  using CxxWrap
  @wrapmodule(joinpath("../../../../build/examples/md-julia/","libjulia_bindings.so"), :define_module_properties)

  function __init__()
    @initcxx
  end
end

module AutoPasM
  using CxxWrap
  @wrapmodule(joinpath("../../../../build/examples/md-julia/","libjulia_bindings.so"), :define_module_autopas)

  function __init__()
    @initcxx
  end
  export updateContainer
end


module ff

include("Configuration.jl")
export
  parseInput

include("ForceCalculation.jl")
export
  updateForces

include("Generator.jl")
export
  generateObject
  generateCubeGrid

include("InputConfiguration.jl")
export
  CubeGridInput
  Thermostat
  InputParameters

include("Simulator.jl")
export
  startSimulation

include("TimeDiscretization.jl")
export
  updatePositions
  updateVelocities

end