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
println("END")