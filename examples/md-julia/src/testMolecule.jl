module autopas
  using CxxWrap
  # @wrapmodule(joinpath("/mnt/c/Users/laura/Documents/BA_SH/autopas_/AutoPas/build/examples/md-julia/","libjulia_bindings.so"))
  @wrapmodule(joinpath("../../../build/examples/md-julia/","libjulia_bindings.so"))

  function __init__()
    @initcxx
  end
end

println("START")

pos = Vector{Float32}([1.1, 2.1, 3.1])
v = Vector{Float64}([1.1, 2.2, 3.3])
molecule = autopas.MoleculeJ{Float32}()
autopas.setPos(molecule, pos)
get_pos = autopas.getPos(molecule)
oldF_set = Vector{Float32}([1.2, 2.3, 3.5])
autopas.setOldF(molecule, oldF_set)

println("type of oldF: " *string(typeof(autopas.getOldF(molecule))))
println("type of pos: " * string(typeof(get_pos)))
println(string(autopas.toString(molecule)))

data_layout = Vector{autopas.DataLayoutOptionValue}([autopas.aos, autopas.soa])
autoPasContainer = autopas.AutoPas{autopas.MoleculeJ{Float64}}()

autopas.setBoxMin(autoPasContainer, [0.0, 0.0, 0.0])
autopas.setBoxMax(autoPasContainer, [7.5, 7.5, 7.5])
autopas.init(autoPasContainer)
println("created container")
autopas.setAllowedDataLayouts(autoPasContainer, data_layout)
it = autopas.IteratorBehavior(autopas.owned)
iter = autopas.begin(autoPasContainer, it)

println("END")