module alib
  using CxxWrap
  # @wrapmodule(joinpath("/mnt/c/Users/laura/Documents/BA_SH/autopas_/AutoPas/build/examples/md-julia/","libjulia_bindings.so"))
  @wrapmodule(joinpath("../../../build/examples/md-julia/","libjulia_bindings.so"))

  function __init__()
    @initcxx
  end
end

# AutoPas/example/md-julia/src

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

# update the velocity of each particle in autoPasContainer
# TODO: where to set iterator behavior and adjust signature (where is deltaT and mass set?)
function updateVelocities(autoPasContainer, deltaT, mass)
  iter = alib.begin(autoPasContainer, alib.owned1)
  while alib.isValid(iter)
    force = alib.getF(alib.deref(iter))
    oldForce = alib.getOldF(alib.deref(iter))
    newV = (oldForce + force) * (deltaT / (2*mass))
    alib.addV(alib.deref(iter), newV)
    alib.inc(iter)

  end
end

# update the positions of each particle in autoPasContainer
# TODO: where to set iterator behavior and adjust signature (where does deltaT, mass and globalForce is set?)
function updatePositions(autoPasContainer, deltaT, mass, globalForce)
  iter = alib.begin(autoPasContainer, alib.owned1)
  while alib.isValid(iter)
    particle = alib.deref(iter)
    velocity = alib.getV(particle)
    force = alib.getF(particle)
    alib.setOldF(particle, force)
    alib.setF(particle, [globalForce, globalForce, globalForce])
    v = velocity * deltaT
    f = force * (deltaT * deltaT / (2*mass))
    vf = v + f
    # TODO: add verlet skin technique and check if particles are too fast
    alib.addPos(particle, vf)
    alib.inc(iter)
  end
end

# update the forces between all particles; iterate over all particles and check the id
function updateForces(autoPasContainer, globalForce, epsilon, sigma)

# calculate pairwise forces: dummy version
  iterOuter = alib.begin(autoPasContainer, alib.owned1)
  while alib.isValid(iterOuter)
    iterInner = alib.begin(autoPasContainer, alib.owned1)
    particleOuter = alib.deref(iterOuter)
    while alib.isValid(iterInner)
      particleInner = alib.deref(iterInner)
      if alib.getID(particleInner) != alib.getID(particleOuter)
        # -24epsilon * 1/||xi-xj||^2 
        # (sigma/||xi-xj||)^6
        # -2*(sigma / |xi-xj||)^12
        # xj-xi
        # diff = broadcast(-, alib.getPos(particleOuter), alib.getPos(particleInner))
        diff = alib.getPos(particleOuter) - alib.getPos(particleInner)  # xi - xj
        pos_diff = (-1) * diff
        # diff = boradcast(*, diff, diff) 
        diff = diff * diff # y:= norm^2 = xik * xjk
        norm = sum(diff) # y0 + y1 + y2
        first = (-24 * epsilon) / norm
        sig6 = (sigma * sigma * sigma * sigma * sigma * sigma)
        norm3 = norm * norm * norm
        second =  sig6 / norm3
        third = -2 * (sig6 * sig6) / (norm3 * norm3)
        scalars = first * second * third
        # force = broadcast(*, particleOuter, constants)
        force = pos_diff * constants
        alib.addF(particleOuter, force)
        # alib.addF(particleInner, broadcast(*, force, -1))
        alib.addF(particleInner, ((-1) * force))
      end
      alib.inc(iterInner)
    end
    alib.inc(iterOuter)
  end

# add global force
  iter = alib.begin(autoPasContainer, alib.owned1)
  while alib.isValid(iter)
    gf = [globalForce, globalForce, globalForce]
    alib.addF(alib.deref(iter), gf)
    alib.inc(iter)
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


h1Hydrogen = alib.()

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
iter = alib.begin(b, alib.owned1)
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

println("start updating velocities")
# container, deltaT, mass
updateVelocities(b, 0.0018, 1.0)
println("end updating velocities")

updatePositions(b, 0.0018, 1.0, 0.0012)
updateForces(b, 0.0018, 1, 1)

#=
for it_ = alib.begin(b), alib.isValid(it_), alib.inc(it_)
  println(alib.toString(alib.deref(it_)))
end
=#

v1 = alib.ContainerOption(alib.directSum)
v2 = alib.ContainerOption(alib.linkedCells)

println("type of v1: " * string(typeof(v1)))
println("type of enum: " * string(typeof(alib.directSum)))

# vec0 = Vector{alib.ContainerOptionValue}([alib.directSum, alib.linkedCells, alib.octree, alib.directSum])
vec0 = Vector{alib.ContainerOption}([v1, v2])

println("type of vector: " * string(typeof(vec0)))
# alib.setAllowedContainers(b, [alib.directSum])
println(typeof(alib.getContainerOption()))
alib.setAllowedContainers(b, vec0)


# alib.setAcquisitionFunction(b, alib.AcquisitionFunctionOption(alib.mean))
#=

ERROR: LoadError: MethodError: no method matching setAllowedContainers(::Main.alib.AutoPasAllocated{Main.alib.MoleculeLJ{Float64}}, ::Vector{Main.alib.ContainerOptionAllocated})
Closest candidates are:
  setAllowedContainers(::Union{CxxWrap.CxxWrapCore.CxxRef{<:Main.alib.AutoPas{Main.alib.MoleculeLJ{Float64}}}, Union{CxxWrap.CxxWrapCore.SmartPointer{T2}, T2} where T2<:Main.alib.AutoPas{Main.alib.MoleculeLJ{Float64}}}, ::Vector{CxxWrap.CxxWrapCore.CxxRef{Main.alib.ContainerOption}}) at ~/.julia/packages/CxxWrap/IdOJa/src/CxxWrap.jl:618
Stacktrace:
 [1] top-level scope
   @ /mnt/c/Users/laura/Documents/BA_SH/AutoPasJulia/AutoPas/examples/md-julia/src/test.jl:240

=#
pi = alib.IteratorBehavior(alib.owned)
iter22 = alib.begin2(b, pi)
res22 = alib.isValid(iter22)
println("is valud: " * string(res22))
println("END")
