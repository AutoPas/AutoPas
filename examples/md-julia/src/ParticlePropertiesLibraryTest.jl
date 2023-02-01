module Properties
  using CxxWrap
  # @wrapmodule(joinpath("/mnt/c/Users/laura/Documents/BA_SH/autopas_/AutoPas/build/examples/md-julia/","libjulia_bindings.so"))
  @wrapmodule(joinpath("../../../build/examples/md-julia/","libjulia_bindings.so"), :define_module_particleLibrary)

  function __init__()
    @initcxx
  end
end
using .Properties

println("START")

p = Properties.ParticlePropertiesLibrary{Float64, Int32}(1.5)
println("created properties object")

Properties.addType(p, 0, 1.2, 1.5, 2.5)
println("added type")

Properties.calculateMixingCoefficients(p)
println("calculated coefficients")

eps = Properties.get24Epsilon(p, convert(Int32, 0));

# md = Properties.getMixingData(p, 0, 0)

println(eps)

println("END")