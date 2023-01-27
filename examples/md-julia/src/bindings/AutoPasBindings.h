#pragma once
#include <iostream>
#include "autopas/AutoPas.h"
/**
 * defining Wrapper for AutoPas class and defining sette
 */

struct WrapAutoPas {
    template<typename Particle>
    void operator()(Particle&& wrapped) {
        // typedef typename T::type WrappedT;
        using WrappedT = typename Particle::type;
        using iterator_t = typename autopas::IteratorTraits<typename WrappedT::Particle_t>::iterator_t;
        wrapped.template constructor<>();
        wrapped.method("init", &WrappedT::init);
        // resizeBox
        // force retune
        wrapped.method("finalize", &WrappedT::finalize);
        wrapped.method("updateContainer", &WrappedT::updateContainer);
        wrapped.method("addParticle", &WrappedT::addParticle);
        // addhaloParticle
        wrapped.method("deleteAllParticles", &WrappedT::deleteAllParticles);
        wrapped.method("begin", static_cast<iterator_t (WrappedT::*)(autopas::options::IteratorBehavior)> (&WrappedT::begin));
        // wrapped.method("deleteParticle", &Wrapped::deleteParticle); may not be wrapped

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

/**
 * defining setters which can not be defined directly
 */

/*
 * setter for ContainerOption
 */
void setAllowedContainers(autopas::AutoPas<autopas::MoleculeLJ<double>>& autoPasContainer, jlcxx::ArrayRef<autopas::options::ContainerOption::Value,1> options) {// jlcxx::ArrayRef<ContainerOptionJulia,1> options) {
    std::set<autopas::options::ContainerOption> tmp;
    for(auto it = options.begin(); it != options.end(); it++) {
        tmp.insert(*it);
    }
    for (auto it = tmp.begin(); it != tmp.end(); it++) {
        std::cout << "options: " << *it << ", ";
    }
    std::cout << "\n";
    autoPasContainer.setAllowedContainers(tmp);
}

/*
 * setter for DataLayout
 */
void setAllowedDataLayouts(autopas::AutoPas<MoleculeJ<double>>& autoPasContainer, jlcxx::ArrayRef<autopas::options::DataLayoutOption::Value,1> options) {
    std::set<autopas::options::DataLayoutOption> tmp;
    for(auto it = options.begin(); it != options.end(); it++) {
        tmp.insert(*it);
    }
    for (auto it = tmp.begin(); it != tmp.end(); it++) {
        std::cout << "options: " << *it << ", ";
    }
    std::cout << "\n";
    autoPasContainer.setAllowedDataLayouts(tmp);
}

/*
 * setter for Newton3Options
 */
void setAllowedNewton3Options(autopas::AutoPas<MoleculeJ<double>>& autoPasContainer, jlcxx::ArrayRef<autopas::options::Newton3Option::Value,1> options) {
    std::set<autopas::options::Newton3Option> tmp;
    for(auto it = options.begin(); it != options.end(); it++) {
        tmp.insert(*it);
    }
    for (auto it = tmp.begin(); it != tmp.end(); it++) {
        std::cout << "options: " << *it << ", ";
    }
    autoPasContainer.setAllowedNewton3Options(tmp);
}

/*
 * setter for TraversalsOption
 */
void setAllowedTraversals(autopas::AutoPas<autopas::MoleculeLJ<double>>& autoPasContainer, jlcxx::ArrayRef<autopas::options::TraversalOption::Value,1> options) {
    std::set<autopas::options::TraversalOption> tmp;
    for(auto it = options.begin(); it != options.end(); it++) {
        tmp.insert(*it);
    }
    for (auto it = tmp.begin(); it != tmp.end(); it++) {
        std::cout << "options: " << *it << ", ";
    }
    std::cout << "\n";
    autoPasContainer.setAllowedTraversals(tmp);
}

/*
 * setter for LoadEstimatorsOption
 */
void setAllowedLoadEstimators(autopas::AutoPas<autopas::MoleculeLJ<double>>& autoPasContainer, jlcxx::ArrayRef<autopas::LoadEstimatorOption::Value,1> options) {
    std::set<autopas::LoadEstimatorOption> tmp;
    for(auto it = options.begin(); it != options.end(); it++) {
        tmp.insert(*it);
    }
    for (auto it = tmp.begin(); it != tmp.end(); it++) {
        std::cout << "options: " << *it << ", ";
    }
    std::cout << "\n";
    autoPasContainer.setAllowedLoadEstimators(tmp);
}

/**
 * setter for boxMin
 */
void setBoxMin(autopas::AutoPas<autopas::MoleculeLJ<double>>& autoPasContainer, jlcxx::ArrayRef<double,1> boxMin) {
    autoPasContainer.setBoxMin({boxMin[0], boxMin[1], boxMin[1]});
}

/**
 * setter for boxMax
 */
void setBoxMax(autopas::AutoPas<autopas::MoleculeLJ<double>>& autoPasContainer, jlcxx::ArrayRef<double,1> boxMax) {
    autoPasContainer.setBoxMax({boxMax[0], boxMax[1], boxMax[1]});
}