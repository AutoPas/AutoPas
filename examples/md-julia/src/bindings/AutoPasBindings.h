#pragma once

struct WrapAutoPas {
    using iterator_t = typename autopas::IteratorTraits<autopas::MoleculeLJ<double>>::iterator_t;
    template<typename Particle>
    void operator()(Particle&& wrapped) {
        // typedef typename T::type WrappedT;
        using WrappedT = typename Particle::type;
        wrapped.template constructor<>();
        wrapped.method("init", &WrappedT::init);
        // resizeBox
        // force retune
        wrapped.method("finalize", &WrappedT::finalize);
        wrapped.method("updateContainer", &WrappedT::updateContainer);
        wrapped.method("addParticle", &WrappedT::addParticle);
        // addhaloParticle
        wrapped.method("deleteAllParticles", &WrappedT::deleteAllParticles);
        wrapped.method("begin2", static_cast<iterator_t (WrappedT::*)(autopas::options::IteratorBehavior)> (&WrappedT::begin));
        // wrapped.method("deleteParticle", &Wrapped::deleteParticle); may not be wrapped
        wrapped.method("printBoxSize", &WrappedT::printBoxSize);

        // setters for AutoPas variables (further setters are defined below)
        wrapped.method("setVerletSkinPerTimestep", &WrappedT::setVerletSkinPerTimestep);
        wrapped.method("setVerletRebuildFrequency", &WrappedT::setVerletRebuildFrequency);
        wrapped.method("setVerletClusterSize", &WrappedT::setVerletClusterSize);
        wrapped.method("setTuningInterval", &WrappedT::setTuningInterval);
        wrapped.method("setMaxEvidence", &WrappedT::setMaxEvidence);
        wrapped.method("setNumSamples", &WrappedT::setNumSamples);
        wrapped.method("setEvidenceFirstPrediction", &WrappedT::setEvidenceFirstPrediction);
        wrapped.method("setRelativeBlacklistRange", &WrappedT::setRelativeBlacklistRange);
        wrapped.method("setMaxTuningPhasesWithoutTest", &WrappedT::setMaxTuningPhasesWithoutTest);
        wrapped.method("setRelativeOptimumRange", &WrappedT::setRelativeOptimumRange);
        wrapped.method("setCutoff", &WrappedT::setCutoff);
        wrapped.method("setOutputSuffix", &WrappedT::setOutputSuffix);
        wrapped.method("setMPITuningWeightForMaxDensity", &WrappedT::setMPITuningWeightForMaxDensity);
        wrapped.method("setMPITuningMaxDifferenceForBucket", &WrappedT::setMPITuningMaxDifferenceForBucket);
        wrapped.method("setAcquisitionFunction", &WrappedT::setAcquisitionFunction);
        wrapped.method("setMPIStrategy", &WrappedT::setMPIStrategy);
        wrapped.method("setTuningStrategyOption", &WrappedT::setTuningStrategyOption);
        wrapped.method("setSelectorStrategy", &WrappedT::setSelectorStrategy);
        wrapped.method("setExtrapolationMethodOption", &WrappedT::setExtrapolationMethodOption);
    }
};

