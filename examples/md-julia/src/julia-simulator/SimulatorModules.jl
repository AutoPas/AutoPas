module Particles
  using CxxWrap
  @wrapmodule(joinpath("../../../../build/examples/md-julia/","libjulia_bindings.so"), :define_module_particles)

  function __init__()
    @initcxx
  end
end

module Options
  using CxxWrap
  @wrapmodule(joinpath("../../../../build/examples/md-julia/","libjulia_bindings.so"), :define_module_options)

  function __init__()
    @initcxx
  end
end

module Iterators
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

end