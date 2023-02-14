# using ..Simulator.Iterators, ..Simulator.AutoPasInterface

function updateForces(autoPasContainer, globalForce, particlePropertiesLibrary)

    tuningIteration = calculatePairwiseForces(autoPasContainer, particlePropertiesLibrary)
    
    iter = AutoPasInterface.begin(autoPasContainer, IteratorBehavior(ownedOrHalo))
    while isValid(iter)
        particle = Simulator.Iterators.:*(iter)
        Particles.addForce(particle, globalForce)
        Simulator.Iterators.:++(iter)
    end
end

# TODO: implement calculatePairwiseForces
function calculatePairwiseForces(autoPasContainer, particlePropertiesLibrary)
    # println("in calculate pairwise")
    return iteratePairwise(autoPasContainer, particlePropertiesLibrary)
end

