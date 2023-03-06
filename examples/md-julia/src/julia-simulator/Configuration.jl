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


function parseInput(inputParameters, comm)

    # create ParticlePropertiesLibrary and generate particles
    particlePropertiesLibrary = Properties.ParticlePropertiesLibrary{Float64, UInt}(inputParameters.cutoff)
    particles = initPplAndParticles(particlePropertiesLibrary, inputParameters.objects)

    # create AutoPas container and initialize all values
    autoPasContainer = initAutoPasContainer(inputParameters)
    Properties.calculateMixingCoefficients(particlePropertiesLibrary)

    # create subdomain
    globalBoxMin = inputParameters.boxMin
    globalBoxMax = inputParameters.boxMax
    
    # globalBoxMin = [0,0,0]
    # globalBoxMax = [10, 10,10]
    # width = calculateMPIWidth(getBoxMin(autoPasContainer), getBoxMax(autoPasContainer), comm)
    width = calculateMPIWidth(globalBoxMin, globalBoxMax, comm)
    # domain = createDomain(getBoxMin(autoPasContainer), getBoxMax(autoPasContainer), width, comm)
    # width = 10
    domain = createDomain(globalBoxMin, globalBoxMax, width, comm, inputParameters.cutoff)

    # add particles to AutoPas container and check if particle is inside local domain
    index = 0
    for particle in particles
        if insideLocalDomain(domain.localBoxMin, domain.localBoxMax, particle)
            addParticle(autoPasContainer, particle)
            index += 1
        end
    end
    # println("added " * string(index) * " particles")
    # rank_id = MPI.Comm_rank(comm)
    println("$domain")

    # println("particles in container: ", getNp(autoPasContainer))

    return autoPasContainer, particlePropertiesLibrary, domain
end

function initPplAndParticles(ppl, objects)
    addedTypes = Dict{Int64, Bool}()
    particles = []
    index = 0
    for particleObject in objects
        
        if !haskey(addedTypes, particleObject.particleType)
            Properties.addType(ppl, particleObject.particleType, particleObject.particleEpsilon, particleObject.particleSigma, particleObject.particleMass)
            addedTypes[particleObject.particleType] = true
        end
        
        # Properties.addType(ppl, particleObject.particleType, particleObject.particleEpsilon, particleObject.particleSigma, particleObject.particleMass)

        generatedParticles, index = generateObject(particleObject, index)
        append!(particles, generatedParticles)
    end
    # println("type of vector: ", typeof(particles), " and length: ", length(particles))
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
    # setVerletRebuildFrequency(apc, inputParameters.verletRebuildFrequency)
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