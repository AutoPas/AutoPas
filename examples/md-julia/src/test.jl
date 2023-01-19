module alib
  using CxxWrap
  @wrapmodule(joinpath("/mnt/c/Users/laura/Documents/BA_SH/AutoPas/build/examples/md-julia/","libjulia_bindings.so"))

  function __init__()
    @initcxx
  end
end

# reason is maybe that some lib is loaded after the oder one?

function generate_dimension(dim)

end


# defining an easy generator for 1D 
function generateParticlesEasy()
  particles = []
  deltaX = 1.12
  xStart = 0
  pos = [0.0, 0.0, 0.0]
  for index in 1:5
    pos[1] = xStart + index * deltaX
    m_ = alib.MoleculeLJ{Float64}()
    alib.setPos(m_, pos)
    push!(particles, m_)
  end
  return particles
end

# generator for a cube with #dimension particles in each dimension 
function generateCube(delta, startPos, dimension)
  particles = []
  for xp in 1:dimension
    xPos = startPos[1] + delta * (xp-1)
    for yp in 1:dimension
      yPos = startPos[2] + delta * (yp-1)
      for zp in 1:dimension
        zPos = startPos[3] + delta * (zp-1)
        m_ = alib.MoleculeLJ{Float64}()
        alib.setPos(m_, [xPos, yPos, zPos])
        push!(particles, m_)
      end
    end
  end
  return particles
end

# print vector of particles (for test purpose)
function printParticles(particles)
  for particle in particles
    println(alib.toString(particle))
  end
end

println("START")

println("START generating particles")
particles = generateCube(1.12, [0.0, 0.0, 0.0], 2)
println("DONE generating particles")

println(typeof(particles))

println("START printing particles")
printParticles(particles)
println("DONE printing particles")

molecule = alib.MoleculeLJ{Float64}()
println("created moleculeLJ instance")
# molecule_1 = alib.MoleculeLJ{Float64}([1.0,1.0,1.0], [2.0,2.0,2.0], convert(Int32, 1), convert(Int32, 1)) # convert(CxxUint, 1), convert(CxxUint, 1))
molecule_1 = alib.MoleculeLJ{Float64}()
alib.setPos(molecule_1, [1.5, 2.5, 3.5])
alib.setV(molecule_1, [2.1, 2.3, 3.3])
println(alib.toString(molecule_1))
# alib.print_particle(molecule) 

pb = alib.ParticleBase{Float64}()
# alib.print_particle(pb)
println("created particlebase")


h1 = alib.Hydrogen()

# x_ = Vector{Float64}([1.0, 1.0, 1.0])
# v_ = Vector{Float64}([5.0, 5.0, 5.0])

x_ = alib.get_vec(1.0, 1.0, 1.0)
v_ = alib.get_vec(5.0, 5.0, 5.0)


h2 = alib.Hydrogen(10, x_, v_)

# a = alib.AutoPas{alib.Hydrogen}()

b = alib.AutoPas{alib.MoleculeLJ{Float64}}()

println("created autopas instance")


alib.setBox(b)
println("set box")

alib.printBoxSize(b)

alib.init(b)
println("init autopas")

alib.getBox(b)
println("get box")

# alib.addParticle(b, molecule)
# println("added particle")

for p in particles
  alib.addParticle(b, p)
end

alib.updateContainer(b)
println("updated container")
# alib.init(b)
# {autopas.MoleculeLJ{Float32}}
# a = alib.AutoPas{alib.MoleculeLJ{Float64}}()
# autopas = autopas.AutoPas{autopas.MoleculeLJ{Float32}}(1)
println("create iterator")
iter = alib.begin(b)
println("done created iterator")
#=
is_valid = alib.isValid(iter)
println("is iterator valid?: " * string(is_valid))
par_ = alib.deref(iter)
println(alib.toString(par_))
iter = alib.inc(iter)
par_ = alib.deref(iter)
println(alib.toString(par_))
=#
while alib.isValid(iter)
  println(alib.toString(alib.deref(iter)))
  alib.inc(iter)
end

#=
for it_ = alib.begin(b), alib.isValid(it_), alib.inc(it_)
  println(alib.toString(alib.deref(it_)))
end
=#

println("END")
