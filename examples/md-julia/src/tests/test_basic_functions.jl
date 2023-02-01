include("../julia-simulator/SimulatorModules.jl")

using .AutoPasM
using .Options
using .Iterators
using .Particles
using .Properties

println("START")

autoPasContianer = AutoPasM.AutoPas{Particles.MoleculeJ{Float64}}()

println("END")