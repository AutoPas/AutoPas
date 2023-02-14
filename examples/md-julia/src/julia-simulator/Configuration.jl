# using ..Properties, ..AutoPasM, ..AIterators
# using ..Simulator.AutoPasInterface, ..Simulator.Iterators


function getNp(apc)
    np0 = 0
    iter = AutoPasInterface.begin(apc, IteratorBehavior(ownedOrHalo))
    while isValid(iter)
        np0 += 1
        Simulator.Iterators.:++(iter)
    end
    return np0
end


function parseInput(inputParameters)

    # create ParticlePropertiesLibrary and generate particles
    particlePropertiesLibrary = Properties.ParticlePropertiesLibrary{Float64, UInt}(inputParameters.cutoff)
    particles = initPplAndParticles(particlePropertiesLibrary, inputParameters.objects)

    # create AutoPas container and initialize all values
    autoPasContainer = initAutoPasContainer(inputParameters)
    Properties.calculateMixingCoefficients(particlePropertiesLibrary)

    # add particles to AutoPas container
    index = 0
    for particle in particles
        addParticle(autoPasContainer, particle)
        index += 1
    end
    println("added " * string(index) * " particles")

    # println("particles in container: ", getNp(autoPasContainer))

    return autoPasContainer, particlePropertiesLibrary
end

function initPplAndParticles(ppl, objects)

    particles = []
    index = 0
    for particleObject in objects
        Properties.addType(ppl, particleObject.particleType, particleObject.particleEpsilon, particleObject.particleSigma, particleObject.particleMass)
        
        generatedParticles, index = generateObject(particleObject, index)
        append!(particles, generatedParticles)
    end
    return particles
end

function initAutoPasContainer(inputParameters)

    apc = AutoPas{Particles.MoleculeJ{Float64}}()

    # add setter with new type
    # setCellS
    setAllowedContainers(apc, inputParameters.container)
    setAllowedDataLayouts(apc, inputParameters.dataLayout)
    setAllowedNewton3Options(apc, inputParameters.newton3)
    setAllowedTraversals(apc, inputParameters.traversal)
    # AutoPasM.setAllowedLoadEstimators(apc, inputparameters._configuration.loadEstimatorOptions.value);
    setBoxMin(apc, inputParameters.boxMin)
    setBoxMax(apc, inputParameters.boxMax)
    setCutoff(apc, inputParameters.cutoff)
    # AutoPasM.setRelativeOptimumRange(_configuration.relativeOptimumRange.value);
    # AutoPasM.setMaxTuningPhasesWithoutTest(_configuration.maxTuningPhasesWithoutTest.value);
    # AutoPasM.setRelativeBlacklistRange(_configuration.relativeBlacklistRange.value);
    # AutoPasM.setEvidenceFirstPrediction(_configuration.evidenceFirstPrediction.value);
    # AutoPasM.setExtrapolationMethodOption(_configuration.extrapolationMethodOption.value);
    # AutoPasM.setNumSamples(_configuration.tuningSamples.value)
    # AutoPasM.setMaxEvidence(_configuration.tuningMaxEvidence.value)
    setSelectorStrategy(apc, inputParameters.selectorStrategy)
    setTuningInterval(apc, inputParameters.tuningInterval)
    setTuningStrategyOption(apc, inputParameters.tuningStrategy)
    setMPIStrategy(apc, inputParameters.mpiStrategy)
    # AutoPasM.setMPITuningMaxDifferenceForBucket(_configuration.MPITuningMaxDifferenceForBucket.value);
    # AutoPasM.setMPITuningWeightForMaxDensity(_configuration.MPITuningWeightForMaxDensity.value);
    # AutoPasM.setVerletClusterSize(apc, inputParameters.verletClusterSize)
    # AutoPasM.setVerletRebuildFrequency(apc, inputParameters.verletRebuildFrequency)
    # AutoPasM.setVerletSkinPerTimestep(apc, inputParameters.verletSkinRadiusPerTimestep)
    # AutoPasM.setAcquisitionFunction(_configuration.acquisitionFunctionOption.value);
    init(apc)
    return apc
end