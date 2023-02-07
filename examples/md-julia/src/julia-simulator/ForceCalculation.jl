function updateForces(autoPasContainer, globalForce)
    # TODO: not implemented yet
    tuningIteration = calculatePairwiseForces(autoPasContainer)

    iter = AutoPasM.begin(autoPasContainer)
    while Iterators.isValid(iter)
        particle = Iterators.:*(iter)
        Particles.addF(particle, gloablForce)
        Iterators.:++(iter)
    end

end

# TODO: implement calculatePairwiseForces
function calculatePairwiseForces(autoPasContainer)
    return false
end

