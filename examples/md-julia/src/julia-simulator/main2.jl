include("./SimulatorModules.jl")
using .Simulator.Particles

function particlePrintTest()
    println("this is some julia test")
    particle = MoleculeJ{Float64}([1.1, 2.2, 3.3], [1.1, 1.1, 1.1], 0, 0)
    println(toString(particle))
end

particlePrintTest()