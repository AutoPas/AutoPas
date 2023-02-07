#include <jlcxx/jlcxx.hpp>
#include <jlcxx/stl.hpp>
#include <type_traits>
#include <jlcxx/functions.hpp>
#include "autopas/molecularDynamics/MoleculeLJ.h"
#include "autopas/molecularDynamics/LJFunctorAVX.h"
#include "TypeDef.h"
#include "autopas/AutoPas.h"
#include "Options.h"
#include "MoleculeJ.h"
#include <iostream>
#include "./bindings/Bindings.h"
/*
#include "Hydrogen.h"
#include <iostream>
#include "IteratorBehaviorJulia.h"
#include "Options.h"
#include "MoleculeJ.h"
#include "./bindings/ParticleBinding.h"
#include "./bindings/AutoPasBindings.h"
#include "./bindings/IteratorBindings.h"
#include "./bindings.AutoPasOptionsBindings.h"
#include "./bindings/ParticlePropertiesLibraryBindings.h"
*/
/*
  numberset _autoPasContainer->setAllowedCellSizeFactors(*_configuration.cellSizeFactors.value);
  DONE set:_autoPasContainer->setAllowedContainers(_configuration.containerOptions.value);
  DONE set _autoPasContainer->setAllowedDataLayouts(_configuration.dataLayoutOptions.value);
  DONE set _autoPasContainer->setAllowedNewton3Options(_configuration.newton3Options.value);
  DONE set _autoPasContainer->setAllowedTraversals(_configuration.traversalOptions.value);
  DONE set _autoPasContainer->setAllowedLoadEstimators(_configuration.loadEstimatorOptions.value);
  DONE array  _autoPasContainer->setBoxMin(_domainDecomposition->getLocalBoxMin());
  DONE array  _autoPasContainer->setBoxMax(_domainDecomposition->getLocalBoxMax());
  DONE O:_autoPasContainer->setCutoff(_configuration.cutoff.value);
  DONE O:_autoPasContainer->setRelativeOptimumRange(_configuration.relativeOptimumRange.value);
  DONE O:_autoPasContainer->setMaxTuningPhasesWithoutTest(_configuration.maxTuningPhasesWithoutTest.value);
  DONE O:_autoPasContainer->setRelativeBlacklistRange(_configuration.relativeBlacklistRange.value);
  DONE O:_autoPasContainer->setEvidenceFirstPrediction(_configuration.evidenceFirstPrediction.value);
  DONE x: _autoPasContainer->setExtrapolationMethodOption(_configuration.extrapolationMethodOption.value);
  DONE O:_autoPasContainer->setNumSamples(_configuration.tuningSamples.value);
  DONE O:_autoPasContainer->setMaxEvidence(_configuration.tuningMaxEvidence.value);
  DONE x: _autoPasContainer->setSelectorStrategy(_configuration.selectorStrategy.value);
  DONE O:_autoPasContainer->setTuningInterval(_configuration.tuningInterval.value);
  DONE x: _autoPasContainer->setTuningStrategyOption(_configuration.tuningStrategyOption.value);
  DONE x: _autoPasContainer->setMPIStrategy(_configuration.mpiStrategyOption.value);
  DONE x: _autoPasContainer->setMPITuningMaxDifferenceForBucket(_configuration.MPITuningMaxDifferenceForBucket.value);
  DONE x: _autoPasContainer->setMPITuningWeightForMaxDensity(_configuration.MPITuningWeightForMaxDensity.value);
  DONE O:_autoPasContainer->setVerletClusterSize(_configuration.verletClusterSize.value);
  DONE O:_autoPasContainer->setVerletRebuildFrequency(_configuration.verletRebuildFrequency.value);
  DONE O:_autoPasContainer->setVerletSkinPerTimestep(_configuration.verletSkinRadiusPerTimestep.value);
  DONE x:_autoPasContainer->setAcquisitionFunction(_configuration.acquisitionFunctionOption.value);
  int rank{};
  autopas::AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &rank);
  DONE O:_autoPasContainer->setOutputSuffix("Rank" + std::to_string(rank) + "_");
  autopas::Logger::get()->set_level(_configuration.logLevel.value);
  DONE  _autoPasContainer->init();
   12/27 total 27 setters
   /1 autopas_mpi_comm_rank
   /1 logger_ set_level
   1/1 init

*/

/**
 * This class helps to wrap the AutoPas type
 */
struct WrapAutoPas {
    template<typename Particle>
    void operator()(Particle&& wrapped) {
        
        using WrappedT = typename Particle::type;
        using iterator_t = typename autopas::IteratorTraits<typename WrappedT::Particle_t>::iterator_t;
        // default constructor 
        wrapped.template constructor<>();
        
        // init method initializing e.g. _tuningStrategy, _autoTuner and _logicHandler
        wrapped.method("init", &WrappedT::init);
        // resizeBox
        // force retune

        // delete particle form container
        wrapped.method("deleteParticle", static_cast<void (WrappedT::*)(typename WrappedT::Particle_t&)> (&WrappedT::deleteParticle));
        
        wrapped.method("finalize", &WrappedT::finalize);

        wrapped.method("updateContainer", &WrappedT::updateContainer);
        
        // add particle to AutoPas container
        wrapped.method("addParticle", &WrappedT::addParticle);
        
        // delete all particles in container 
        wrapped.method("deleteAllParticles", &WrappedT::deleteAllParticles);
        
        // return iterator to first element of container
        wrapped.method("begin", static_cast<iterator_t (WrappedT::*)(autopas::options::IteratorBehavior)> (&WrappedT::begin));

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

        // getters of AutoPas variables
        wrapped.method("getCutoff", &WrappedT::getCutoff);    
    }
};

/**
 * wrapper method for iteratePairwise
 */
/*
template<class Functor>
bool iteratePairwise(autopas::AutoPas<MoleculeJ<double>>& autoPasContainer, Functor *f) {
    return autoPasContainer.iteratePairwise(f);
}
*/
/**
 * wrapper method for iteratePairwise with constant functor LJFunctorAVX
 
bool iteratePairwise(autopas::AutoPas<MoleculeJ<double>>& autoPasContainer, autopas::LJFunctorAVX<MoleculeJ<double>>* functor) {
    return autoPasContainer.iteratePairwise(functor);
}
*/

/**
 * create functor LJFunctorAVX and return it
 */
/*
autopas::LJFunctorAVX<MoleculeJ<double>> getLJFunctorAVX(autopas::AutoPas<MoleculeJ<double>>& autoPasContainer, ParticlePropertiesLibrary<double, size_t>& particlePropertiesLibrary) {
    autopas::LJFunctorAVX<MoleculeJ<double>, true, true> functor{autoPasContainer.getCutoff(), particlePropertiesLibrary};
    return functor;
}
*/

/**
 * wrapper methods for setters which can not be defined directly
 */

/*
 * setter for ContainerOption
 */
void setAllowedContainers(autopas::AutoPas<MoleculeJ<double>>& autoPasContainer, jlcxx::ArrayRef<autopas::options::ContainerOption::Value,1> options) {// jlcxx::ArrayRef<ContainerOptionJulia,1> options) {
    std::set<autopas::options::ContainerOption> tmp;
    for(auto it = options.begin(); it != options.end(); it++) {
        tmp.insert(*it);
    }
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
    autoPasContainer.setAllowedNewton3Options(tmp);
}

/*
 * setter for TraversalsOption
 */
void setAllowedTraversals(autopas::AutoPas<MoleculeJ<double>>& autoPasContainer, jlcxx::ArrayRef<autopas::options::TraversalOption::Value,1> options) {
    std::set<autopas::options::TraversalOption> tmp;
    for(auto it = options.begin(); it != options.end(); it++) {
        tmp.insert(*it);
    }
    autoPasContainer.setAllowedTraversals(tmp);
}

/*
 * setter for LoadEstimatorsOption
 */
void setAllowedLoadEstimators(autopas::AutoPas<MoleculeJ<double>>& autoPasContainer, jlcxx::ArrayRef<autopas::LoadEstimatorOption::Value,1> options) {
    std::set<autopas::LoadEstimatorOption> tmp;
    for(auto it = options.begin(); it != options.end(); it++) {
        tmp.insert(*it);
    }
    autoPasContainer.setAllowedLoadEstimators(tmp);
}

/**
 * setter for boxMin
 */
void setBoxMin(autopas::AutoPas<MoleculeJ<double>>& autoPasContainer, jlcxx::ArrayRef<double,1> boxMin) {
    autoPasContainer.setBoxMin({boxMin[0], boxMin[1], boxMin[1]});
}

/**
 * setter for boxMax
 */
void setBoxMax(autopas::AutoPas<MoleculeJ<double>>& autoPasContainer, jlcxx::ArrayRef<double,1> boxMax) {
    autoPasContainer.setBoxMax({boxMax[0], boxMax[1], boxMax[1]});
}

/**
 * create Julia Module with AutoPas type
 */
JLCXX_MODULE define_module_autopas(jlcxx::Module& mod)
{
    using jlcxx::Parametric;
    using jlcxx::TypeVar;

    /**
     * add AutoPas tyoe to Julia with template parameter
     * autopas::MoleculeLJ<double> and MoleculeJ<double>
     */
    mod.add_type<Parametric<TypeVar<1>>>("AutoPas")
            .apply<autopas::AutoPas<autopas::MoleculeLJ<double>>, autopas::AutoPas<MoleculeJ<double>>>(WrapAutoPas());

    /**
     * add wrapper method for iteratePairwise of AutoPas
     */
    // mod.method("iteratePairwise", &iteratePairwise);

    /**
     * add wrapper method for getLJFunctorAVX to Julia
     */
    // mod.method("getLJFunctorAVX", &getLJFunctorAVX);
    
    /*
     * add wrapper methods for setters of AutoPas attributes which cannot be wrapped directly
     */
    
    /*
     * add setAllowedContainers method to Julia
     */
    mod.method("setAllowedContainers", &setAllowedContainers);

    /*
     * add setAllowedDataLayouts method to Julia
     */
    mod.method("setAllowedDataLayouts", setAllowedDataLayouts);

    /*
     * add setAllowedNewton3Options method to Julia
     */
    mod.method("setAllowedNewton3Options", setAllowedNewton3Options);

    /*
     * add setAllowedTraversals method to Julia
     */
    mod.method("setAllowedTraversals", setAllowedTraversals);

    /*
     * add setAllowedLoadEstimators method to Julia
     */
    mod.method("setAllowedLoadEstimators", setAllowedLoadEstimators);

    /**
     * add setBoxMin method to Julia
     */
    mod.method("setBoxMin", &setBoxMin);

    /**
     * add setBoxMax method to Julia
     */
    mod.method("setBoxMax", &setBoxMax);
}

// cmake -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DAUTOPAS_BUILD_TESTS=OFF ..