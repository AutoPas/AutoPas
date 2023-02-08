using ..Properties, ..AutoPasM, ..Iterators

function parseInput(inputParameters)

    # create ParticlePropertiesLibrary and generate particles
    particlePropertiesLibrary = Properties.ParticlePropertiesLibrary{Float64, UInt}(inputParameters.cutoff)
    particles = initPplAndParticles(particlePropertiesLibrary, inputParameters.objects)

    # create AutoPas container and initialize all values
    autoPasContainer = initAutoPasContainer(inputParameters)
    Properties.calculateMixingCoefficients(particlePropertiesLibrary)

    # add particles to AutoPas container
    for particle in particles
        AutoPasM.addParticle(autoPasContainer, particle)
    end

    iter = AutoPasM.begin(autoPasContainer, Options.IteratorBehavior(Options.owned))
    while Iterators.isValid(iter)
        println(Particles.toString(Iterators.:*(iter)))
        Iterators.:++(iter)
    end
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

    apc = AutoPasM.AutoPas{Particles.MoleculeJ{Float64}}()

    # add setter with new type
    # setCellS
    AutoPasM.setAllowedContainers(apc, inputParameters.container)
    AutoPasM.setAllowedDataLayouts(apc, inputParameters.dataLayout)
    AutoPasM.setAllowedNewton3Options(apc, inputParameters.newton3)
    AutoPasM.setAllowedTraversals(apc, inputParameters.traversal)
    # AutoPasM.setAllowedLoadEstimators(apc, inputparameters._configuration.loadEstimatorOptions.value);
    AutoPasM.setBoxMin(apc, inputParameters.boxMin)
    AutoPasM.setBoxMax(apc, inputParameters.boxMax)
    AutoPasM.setCutoff(apc, inputParameters.cutoff)
    # AutoPasM.setRelativeOptimumRange(_configuration.relativeOptimumRange.value);
    # AutoPasM.setMaxTuningPhasesWithoutTest(_configuration.maxTuningPhasesWithoutTest.value);
    # AutoPasM.setRelativeBlacklistRange(_configuration.relativeBlacklistRange.value);
    # AutoPasM.setEvidenceFirstPrediction(_configuration.evidenceFirstPrediction.value);
    # AutoPasM.setExtrapolationMethodOption(_configuration.extrapolationMethodOption.value);
    # AutoPasM.setNumSamples(_configuration.tuningSamples.value)
    # AutoPasM.setMaxEvidence(_configuration.tuningMaxEvidence.value)
    AutoPasM.setSelectorStrategy(apc, inputParameters.selectorStrategy)
    AutoPasM.setTuningInterval(apc, inputParameters.tuningInterval)
    AutoPasM.setTuningStrategyOption(apc, inputParameters.tuningStrategy)
    AutoPasM.setMPIStrategy(apc, inputParameters.mpiStrategy)
    # AutoPasM.setMPITuningMaxDifferenceForBucket(_configuration.MPITuningMaxDifferenceForBucket.value);
    # AutoPasM.setMPITuningWeightForMaxDensity(_configuration.MPITuningWeightForMaxDensity.value);
    AutoPasM.setVerletClusterSize(apc, inputParameters.verletClusterSize)
    AutoPasM.setVerletRebuildFrequency(apc, inputParameters.verletRebuildFrequency)
    AutoPasM.setVerletSkinPerTimestep(apc, inputParameters.verletSkinRadiusPerTimestep)
    # AutoPasM.setAcquisitionFunction(_configuration.acquisitionFunctionOption.value);
    AutoPasM.init(apc)
    return apc
end