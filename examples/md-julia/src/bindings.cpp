#include <jlcxx/jlcxx.hpp>
#include <jlcxx/stl.hpp>
#include <type_traits>
#include <jlcxx/functions.hpp>
#include "autopas/molecularDynamics/MoleculeLJ.h"
#include "autopas/molecularDynamics/LJFunctorAVX.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
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
    template<typename T>
    void operator()(T&& autoPas) {
        
        using AutoPasType = typename T::type;
        using iterator_t = typename autopas::IteratorTraits<typename AutoPasType::Particle_t>::iterator_t;
        // default constructor 
        autoPas.template constructor<>();
        
        // init method initializing e.g. _tuningStrategy, _autoTuner and _logicHandler
        autoPas.method("init", &AutoPasType::init);
        // resizeBox
        // force retune

        // delete particle form container
        autoPas.method("deleteParticle", static_cast<void (AutoPasType::*)(typename AutoPasType::Particle_t&)> (&AutoPasType::deleteParticle));
        
        autoPas.method("finalize", &AutoPasType::finalize);

        autoPas.method("updateContainer", &AutoPasType::updateContainer);
        
        // add particle to AutoPas container
        autoPas.method("addParticle", &AutoPasType::addParticle);
        
        // delete all particles in container 
        autoPas.method("deleteAllParticles", &AutoPasType::deleteAllParticles);
        
        // return iterator to first element of container
        autoPas.method("begin", static_cast<iterator_t (AutoPasType::*)(autopas::options::IteratorBehavior)> (&AutoPasType::begin));

        // retrun number of particles in the container
        autoPas.method("getNumberOfParticles", &AutoPasType::getNumberOfParticles);

        // add a halo particle to the container
        autoPas.method("addHaloParticle", &AutoPasType::addHaloParticle);

        // setters for AutoPas variables (further setters are defined below)
        autoPas.method("setVerletSkinPerTimestep", &AutoPasType::setVerletSkinPerTimestep);
        autoPas.method("setVerletRebuildFrequency", &AutoPasType::setVerletRebuildFrequency);
        autoPas.method("setVerletClusterSize", &AutoPasType::setVerletClusterSize);
        autoPas.method("setTuningInterval", &AutoPasType::setTuningInterval);
        autoPas.method("setMaxEvidence", &AutoPasType::setMaxEvidence);
        autoPas.method("setNumSamples", &AutoPasType::setNumSamples);
        autoPas.method("setEvidenceFirstPrediction", &AutoPasType::setEvidenceFirstPrediction);
        autoPas.method("setRelativeBlacklistRange", &AutoPasType::setRelativeBlacklistRange);
        autoPas.method("setMaxTuningPhasesWithoutTest", &AutoPasType::setMaxTuningPhasesWithoutTest);
        autoPas.method("setRelativeOptimumRange", &AutoPasType::setRelativeOptimumRange);
        autoPas.method("setCutoff", &AutoPasType::setCutoff);
        autoPas.method("setOutputSuffix", &AutoPasType::setOutputSuffix);
        autoPas.method("setMPITuningWeightForMaxDensity", &AutoPasType::setMPITuningWeightForMaxDensity);
        autoPas.method("setMPITuningMaxDifferenceForBucket", &AutoPasType::setMPITuningMaxDifferenceForBucket);
        autoPas.method("setAcquisitionFunction", &AutoPasType::setAcquisitionFunction);
        autoPas.method("setMPIStrategy", &AutoPasType::setMPIStrategy);
        autoPas.method("setTuningStrategyOption", &AutoPasType::setTuningStrategyOption);
        autoPas.method("setSelectorStrategy", &AutoPasType::setSelectorStrategy);
        autoPas.method("setExtrapolationMethodOption", &AutoPasType::setExtrapolationMethodOption);

        // getters of AutoPas variables
        autoPas.method("getCutoff", &AutoPasType::getCutoff);
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
 */
bool iteratePairwise(autopas::AutoPas<MoleculeJ<double>>& autoPasContainer, ParticlePropertiesLibrary<> particlePropertiesLibrary) {
    /*std::cout << "in c++ function: ";
    for (auto it = autoPasContainer.begin(); it.isValid(); ++it) {
        std::cout << "[" << (*it).getR()[0] << ", " << (*it).getR()[1] << ", " << (*it).getR()[2] << "], ";
        std::cout << "[" << (*it).getV()[0] << ", " << (*it).getV()[1] << ", " << (*it).getV()[2] << "], ";
        std::cout << "[" << (*it).getF()[0] << ", " << (*it).getF()[1] << ", " << (*it).getF()[2] << "]\n";
        // std::cout << (*it).toString() << " ";
    }
    std::cout << std::endl;
    */
    autopas::LJFunctorAVX<MoleculeJ<double>, true, true> functor{autoPasContainer.getCutoff(), particlePropertiesLibrary};
    autoPasContainer.iteratePairwise(&functor);
    /*
    for (auto it = autoPasContainer.begin(); it.isValid(); ++it) {
        std::cout << "[" << (*it).getR()[0] << ", " << (*it).getR()[1] << ", " << (*it).getR()[2] << "], ";
        std::cout << "[" << (*it).getV()[0] << ", " << (*it).getV()[1] << ", " << (*it).getV()[2] << "], ";
        std::cout << "[" << (*it).getF()[0] << ", " << (*it).getF()[1] << ", " << (*it).getF()[2] << "]\n";
        // std::cout << (*it).toString() << " ";
    }
    std::cout << std::endl;
    */
    return false;
}

/**
 * wrapper method for region iterator
 */
autopas::IteratorTraits<MoleculeJ<double>>::iterator_t regionIterator(autopas::AutoPas<MoleculeJ<double>>& autoPasContainer, jlcxx::ArrayRef<double, 1> lowerCorner, jlcxx::ArrayRef<double,1> upperCorner, autopas::options::IteratorBehavior behavior) {
    std::array<double,3> lC{lowerCorner[0], lowerCorner[1], lowerCorner[2]};
    std::array<double,3> uC{upperCorner[0], upperCorner[1], upperCorner[2]};
    return autoPasContainer.getRegionIterator(lC, uC, behavior);
}

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

jlcxx::ArrayRef<double,1> getBoxMin(autopas::AutoPas<MoleculeJ<double>>& autoPasContainer) {
    // std::array<double,3> tmp = autoPasContainer.getBoxMin();
    return {autoPasContainer.getBoxMin().data(), autoPasContainer.getBoxMin().size()};
}

jlcxx::ArrayRef<double,1> getBoxMax(autopas::AutoPas<MoleculeJ<double>>& autoPasContainer) {
    // std::array<double,3> tmp = autoPasContainer.getBoxMax(); 
    return {autoPasContainer.getBoxMax().data(), autoPasContainer.getBoxMax().size()};
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
    mod.method("iteratePairwise", &iteratePairwise);

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

    /**
     * add getBoxMin mthod to Julia
     */
    mod.method("getBoxMin", &getBoxMin);
    
    /**
     * add getBoxMax mthod to Julia
     */
    mod.method("getBoxMax", &getBoxMax);

    /**
     * add regionIterator to Julia
     */
    mod.method("regionIterator", &regionIterator);
}

// cmake -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DAUTOPAS_BUILD_TESTS=OFF -DAUTOPAS_OPENMP=ON ..
// cmake -DMD_FLEXIBLE_USE_MPI=ON -DMD_FLEXIBLE_FUNCTOR_AUTOVEC=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DAUTOPAS_BUILD_TESTS=OFF -DAUTOPAS_OPENMP=ON ..
// cmake -DMD_FLEXIBLE_USE_MPI=OFF -DMD_FLEXIBLE_FUNCTOR_AUTOVEC=OFF -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DAUTOPAS_BUILD_TESTS=OFF -DAUTOPAS_OPENMP=OFF ..