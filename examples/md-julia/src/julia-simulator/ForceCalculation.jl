function updateForces(autoPasContainer, globalForce, particlePropertiesLibrary)
    # TODO: not implemented yet
    tuningIteration = calculatePairwiseForces(autoPasContainer, particlePropertiesLibrary)
    #=
    iter = AutoPasM.begin(autoPasContainer, Options.IteratorBehavior(Options.owned))
    while AIterators.isValid(iter)
        particle = AIterators.:*(iter)
        Particles.addForce(particle, globalForce)
        AIterators.:++(iter)
    end
    =#

end

# TODO: implement calculatePairwiseForces
function calculatePairwiseForces(autoPasContainer, particlePropertiesLibrary)
    println("in calculate pairwise")
    AutoPasM.iteratePairwise(autoPasContainer, particlePropertiesLibrary)
    return false
end

