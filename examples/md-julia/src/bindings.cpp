#include <jlcxx/jlcxx.hpp>
#include <jlcxx/stl.hpp>
#include <type_traits>
#include "jlcxx/functions.hpp"
#include "autopas/molecularDynamics/MoleculeLJ.h"
// #include "autopas/AutoPas.h"
#include "TypeDef.h"
// #include "Simulation.h"
#include "autopas/AutoPas.h"
#include "Hydrogen.h"
#include <iostream>
#include "IteratorBehaviorJulia.h"
#include "Options.h"
#include "MoleculeJ.h"
#include "./bindings/MoleculeBinding.h"
#include "./bindings/AutoPasBindings.h"


// extern template class autopas::AutoPas<ParticleType>;
// template class autopas::AutoPas<ParticleType>;

struct WrapMoleculeLJ {
    template<typename floatType>
    void operator()(floatType&& wrapped) {
        // typedef typename T::type WrappedT;
        using WrappedT = typename floatType::type;

        // constructors for MoleculeLJ
        wrapped.template constructor<>();
        wrapped.template constructor<jlcxx::ArrayRef<double,1>, jlcxx::ArrayRef<double,1>, int, int>();

        // setters of MoleculeLJ attributes (adjusted for the Julia Wrapper)
        wrapped.method("setPos", &WrappedT::setPos);
        wrapped.method("setV", &WrappedT::setV);
        wrapped.method("setF", &WrappedT::setF);
        wrapped.method("setOldF", static_cast<void (WrappedT::*)(jlcxx::ArrayRef<double,1>)> (&WrappedT::setOldF));

        // getters of MoleculeLJ attributes (adjusted for the Julia Wrapper)
        wrapped.method("getPos", &WrappedT::getPos);
        wrapped.method("getV", &WrappedT::getV);
        wrapped.method("getF", &WrappedT::getF);
        wrapped.method("getOldF", static_cast<jlcxx::ArrayRef<double,1> (WrappedT::*)()> (&WrappedT::getOldF));
//         //   t0.method("set_x_pos", static_cast<void (WrapClass::*)(int) >(&WrapClass::set_x_pos));

        // add methods of MoleculeLJ attributes (adjusted for the Julia Wrapper)
        wrapped.method("addPos", &WrappedT::addPos);
        wrapped.method("addV", &WrappedT::addV);
        wrapped.method("addF", &WrappedT::addF);

        // 
        wrapped.method("getID", &WrappedT::getID);
        wrapped.method("toString", &WrappedT::toString);
        // wrapped.method("print_particle", &WrappedT::print_particle);
    }
};

struct WrapParticleBase {
    template<typename T>
    void operator()(T&& wrapped) {
        using WrapedT = typename T::type;
        wrapped.template constructor<>();
        // wrapped.method("print_particle", &WrapedT::print_particle);
    }
};

struct WrapIteratorInterface {
    template<typename T>
    void operator()(T&& wrapped) {
        using WrapedT = typename T::type;
        // wrapped.template constructor<>();
        
    }
};

/*
struct WrapIteratorInterfaceImpl {
    template<typename T>
    void operator()(T&& wrapped) {
        using WrappedT = typename T::type;
        // wrapped.template constructor<>();
    }

};
*/

struct WrapIteratorWrapper {
    template<typename T>
    void operator()(T&& wrapped) {
        using WrappedT = typename T::type;
        // wrapped.template constructor<>();
        // wrapped.method("setIterator", &WrappedT::setInterator);
        wrapped.method("isValid", &WrappedT::isValid);
        wrapped.method("inc", &WrappedT::operator++);
        wrapped.method("deref", &WrappedT::operator*);
        // TODO: wrap operator
    }
};

namespace jlcxx
{
    using jlcxx::Parametric;
    // specifying BuildParameterList
    
    
    template<typename Particle, bool modifiable>
    struct BuildParameterList<autopas::ParticleIteratorInterface<Particle, modifiable>> {
        typedef ParameterList<Particle, bool> type;
    };

    /*
    template<typename Particle, bool modifiable>
    struct BuildParameterList<autopas::internal::ParticleIteratorInterfaceImpl<Particle, modifiable>> {
        typedef ParameterList<Particle, bool> type;
    };
    */

    template<typename Particle, bool modifiable>
    struct BuildParameterList<autopas::ParticleIteratorWrapper<Particle, modifiable>> {
        typedef ParameterList<Particle, bool> type;
    };
    
    // adding inheritance for types
    template<typename Particle, bool modifiable> struct SuperType<autopas::ParticleIteratorWrapper<Particle, modifiable>> {typedef autopas::ParticleIteratorInterface<Particle, modifiable> type;};

    // adding templates for types
    template<typename floatType> struct IsMirroredType<autopas::MoleculeLJ<floatType>> : std::false_type { };
    template<typename floatType> struct IsMirroredType<MoleculeJ<floatType>> : std::false_type { };
    template<typename Particle> struct IsMirroredType<autopas::AutoPas<Particle>> : std::false_type { };
    template<typename T, typename V> struct IsMirroredType<autopas::ParticleBase<T, V>> : std::false_type { };
    template<typename Particle, bool modifiable> struct IsMirroredType<autopas::ParticleIteratorInterface<Particle, modifiable>> : std::false_type { };
    template<typename Particle, bool modifiable> struct IsMirroredType<autopas::ParticleIteratorWrapper<Particle, modifiable>> : std::false_type { };
}

/*
  numberset _autoPasContainer->setAllowedCellSizeFactors(*_configuration.cellSizeFactors.value);
  DONE set:_autoPasContainer->setAllowedContainers(_configuration.containerOptions.value);
  set _autoPasContainer->setAllowedDataLayouts(_configuration.dataLayoutOptions.value);
  set _autoPasContainer->setAllowedNewton3Options(_configuration.newton3Options.value);
  set _autoPasContainer->setAllowedTraversals(_configuration.traversalOptions.value);
  set _autoPasContainer->setAllowedLoadEstimators(_configuration.loadEstimatorOptions.value);
  array  _autoPasContainer->setBoxMin(_domainDecomposition->getLocalBoxMin());
  array  _autoPasContainer->setBoxMax(_domainDecomposition->getLocalBoxMax());
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
   _autoPasContainer->init();
   12/27 total 27 setters
   /1 autopas_mpi_comm_rank
   /1 logger_ set_level
   1/1 init

*/

/*
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
        // .method("open", static_cast<void (IO::LCReader::*)(const std::string&)>(&IO::LCReader::open));
        //  t0.method("get_x_pos", static_cast<int (WrapClass::*)() >(&WrapClass::get_x_pos));
        //   t0.method("set_x_pos", static_cast<void (WrapClass::*)(int) >(&WrapClass::set_x_pos));
        // wrapped.method("begin", static_cast<iterator_t (WrappedT::*)(autopas::IteratorBehavior) > (&WrappedT::begin));
        // wrapped.method("begin", static_cast<autopas::ParticleIteratorWrapper<Particle, bool> (WrappedT::*)() > (&WrappedT::begin));
        // wrapped.method("begin", &WrappedT::begin);
        wrapped.method("begin2", static_cast<iterator_t (WrappedT::*)(autopas::options::IteratorBehavior)> (&WrappedT::begin));
        // wrapped.method("begin", static_cast<autopas::ParticleIteratorWrapper<autopas::MoleculeLJ<double>, true> (WrappedT::*)() >(&WrappedT::begin));
        // wrapped.method("deleteParticle", &Wrapped::deleteParticle); may not be wrapped
        wrapped.method("printBoxSize", &WrappedT::printBoxSize);

        // setters for AutoPas variables (direct)
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
*/

/*
 * further setters of the AutoPas class which can not be wraped directly
 */

/*
 * setter for ContainerOption
 */
// void setAllowedContainers(autopas::AutoPas<autopas::MoleculeLJ<double>>& autoPasContainer, ContainerOptionJulia op) {// jlcxx::ArrayRef<ContainerOptionJulia,1> option) {
    /*std::set<autopas::options::ContainerOption> tmp;
    for(auto it = option.begin(); it != option.end(); it++) {
        ContainerOptionJulia c = *it;
        // autopas::options::ContainerOption op = it;
        // autopas::options::ContainerOption x = static_cast<autopas::options::ContainerOption>(c);
        tmp.insert(static_cast<autopas::options::ContainerOption>(1));
    }
    */
   // std::set<autopas::options::ContainerOption> tmp;
    // autopas::options::ContainerOption oo = static_cast<int>(op);
    // tmp.insert(static_cast<autopas::options::ContainerOption>(op));
    // autoPasContainer.setAllowedContainers(tmp);
// }

/*
 * setter for DataLayout
 */
 /*
void setAllowedDataLayouts(autopas::AutoPas<autopas::MoleculeLJ<double>>& autoPasContainer, jlcxx::ArrayRef<DataLayoutOptionJulia,1> option) {
    autoPasContainer.setAllowedDataLayouts({option.begin(), option.end()});
}
*/
/*
 * setter for Newton3Options
 */
 /*
void setAllowedNewton3Options(autopas::AutoPas<autopas::MoleculeLJ<double>>& autoPasContainer, jlcxx::ArrayRef<Newton3OptionJulia,1> option) {
    autoPasContainer.setAllowedNewton3Options({option.begin(), option.end()});
}
*/

/*
 * setter for TraversalsOption
 */
/*
void setAllowedTraversals(autopas::AutoPas<autopas::MoleculeLJ<double>>& autoPasContainer, jlcxx::ArrayRef<TraversalOptionJulia,1> option) {
    autoPasContainer.setAllowedTraversals({option.begin(), option.end()});
}
*/

/*
 * setter for LoadEstimatorsOption
 */
/*
void setAllowedLoadEstimators(autopas::AutoPas<autopas::MoleculeLJ<double>>& autoPasContainer, jlcxx::ArrayRef<LoadEstimatorOptionJulia,1> option) {
    autoPasContainer.setAllowedLoadEstimators({option.begin(), option.end()});
}
*/

/*
struct WrapOptions {
    template<typename actualOption>
    void operator()(actualOption&& wrapped) {
        using WrappedT = typename actualOption::type;
    }
};
*/

// void setAllowedContainers(autopas::AutoPas<autopas::MoleculeLJ<double>>& autoPasContainer, jlcxx::ArrayRef<autopas::options::ContainerOption,1> options) {
void setAllowedContainers(autopas::AutoPas<autopas::MoleculeLJ<double>>& autoPasContainer, jlcxx::ArrayRef<autopas::options::ContainerOption,1> option) {// jlcxx::ArrayRef<autopas::options::ContainerOption,1> option) {
    /*std::set<autopas::options::ContainerOption> tmp;
    for(auto it = options.begin(); it != options.end(); it++) {
        tmp.insert(*it);
    }

    for(auto it = tmp.begin(); it != tmp.end(); it++) {
        std::cout << "option: " << *it << ", ";
    }
    std::cout << std::endl;
    autoPasContainer.setAllowedContainers(tmp);
    */
    /*
    for(auto it = option.begin(); it != option.end(); it++) {
        std::cout << "option: " << *it << ", ";
    }
    */
    std::cout << "in set allowed container\n";
}

autopas::options::ContainerOption getContainerOption() {
    autopas::options::ContainerOption option{autopas::options::ContainerOption::directSum};
    return option;
}

std::vector<double> get_vec(double x1, double x2, double x3) {
    return {x1, x2, x3};
}


void setBox(autopas::AutoPas<autopas::MoleculeLJ<double>>& ap) {
    const std::array<double, 3> min_{0,0,0};
    const std::array<double, 3> max_{7.5,7.5,7.5};
    
    ap.setBoxMin(min_);
    ap.setBoxMax(max_);
    /*
    auto m1 = ap.getBoxMin();
    auto m2 = ap.getBoxMax();
    std::cout << "min: " << m1.at(0) << ", " << m1.at(1) << ", " << m1.at(2) << "\n";
    std::cout << "min: " << m2.at(0) << ", " << m2.at(1) << ", " << m2.at(2) << "\n";
    */
}

void getBox(autopas::AutoPas<autopas::MoleculeLJ<double>>& ap) {
    auto m1 = ap.getBoxMin();
    auto m2 = ap.getBoxMax();
    std::cout << "get_min: " << m1.at(0) << ", " << m1.at(1) << ", " << m1.at(2) << "\n";
    std::cout << "get_min: " << m2.at(0) << ", " << m2.at(1) << ", " << m2.at(2) << "\n";
}

autopas::ParticleIteratorWrapper<autopas::MoleculeLJ<double>, true> begin_new(autopas::AutoPas<autopas::MoleculeLJ<double>>& autopasContianer, IteratorBehaviorJulia iteratorBehaviorJulia) {
    autopas::options::IteratorBehavior b = iteratorBehaviorJulia;
    std::cout << "in begin_new iterator function\n";
    return autopasContianer.begin(b);
    // return autopasContianer.begin();
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    using jlcxx::Parametric;
    using jlcxx::TypeVar;

    mod.add_bits<autopas::options::AcquisitionFunctionOption::Value>("AcquisitionFunctionOptionValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("upperConfidenceBound", autopas::options::AcquisitionFunctionOption::Value::upperConfidenceBound);
    mod.set_const("mean", autopas::options::AcquisitionFunctionOption::Value::mean);
    mod.set_const("variance", autopas::options::AcquisitionFunctionOption::Value::variance);
    mod.set_const("probabilityOfImprovement", autopas::options::AcquisitionFunctionOption::Value::probabilityOfImprovement);
    mod.set_const("expectedImprovement", autopas::options::AcquisitionFunctionOption::Value::expectedImprovement);

    mod.add_type<autopas::options::AcquisitionFunctionOption>("AcquisitionFunctionOption")
            .constructor<autopas::options::AcquisitionFunctionOption::Value>();

    mod.add_bits<autopas::options::MPIStrategyOption::Value>("MPIStrategyOptionValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("noMPI", autopas::options::MPIStrategyOption::Value::noMPI);
    mod.set_const("divideAndConquer", autopas::options::MPIStrategyOption::Value::divideAndConquer);

    mod.add_type<autopas::options::MPIStrategyOption>("MPIStrategyOption")
            .constructor<autopas::options::MPIStrategyOption::Value>();

    mod.add_bits<autopas::options::TuningStrategyOption::Value>("TuningStrategyOptionValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("randomSearch", autopas::options::TuningStrategyOption::Value::randomSearch);
    mod.set_const("fullSearch", autopas::options::TuningStrategyOption::Value::fullSearch);
    mod.set_const("bayesianSearch", autopas::options::TuningStrategyOption::Value::bayesianSearch);
    mod.set_const("bayesianClusterSearch", autopas::options::TuningStrategyOption::Value::bayesianClusterSearch);
    mod.set_const("activeHarmony", autopas::options::TuningStrategyOption::Value::activeHarmony);
    mod.set_const("predictiveTuning", autopas::options::TuningStrategyOption::Value::predictiveTuning);

    mod.add_type<autopas::options::TuningStrategyOption>("TuningStrategyOption")
            .constructor<autopas::options::TuningStrategyOption::Value>();

    mod.add_bits<autopas::options::SelectorStrategyOption::Value>("SelectorStrategyOptionValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("fastestAbs", autopas::options::SelectorStrategyOption::Value::fastestAbs);
    mod.set_const("fastestMean", autopas::options::SelectorStrategyOption::Value::fastestMean);
    mod.set_const("fastestMedian", autopas::options::SelectorStrategyOption::Value::fastestMedian);

    mod.add_type<autopas::options::SelectorStrategyOption>("SelectorStrategyOption")
            .constructor<autopas::options::SelectorStrategyOption::Value>();

    mod.add_bits<autopas::options::ExtrapolationMethodOption::Value>("ExtrapolationMethodOptionValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("linePrediction", autopas::options::ExtrapolationMethodOption::Value::linePrediction);
    mod.set_const("linearRegression", autopas::options::ExtrapolationMethodOption::Value::linearRegression);
    mod.set_const("newton", autopas::options::ExtrapolationMethodOption::Value::newton);

    mod.add_type<autopas::options::ExtrapolationMethodOption>("ExtrapolationMethodOption")
            .constructor<autopas::options::ExtrapolationMethodOption::Value>();

    mod.add_bits<autopas::options::IteratorBehavior::Value>("IteratorBehaviorValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("owned", autopas::options::IteratorBehavior::Value::owned);
    mod.set_const("halo", autopas::options::IteratorBehavior::Value::halo);
    mod.set_const("ownedOrHalo", autopas::options::IteratorBehavior::Value::ownedOrHalo);
    mod.set_const("dummy", autopas::options::IteratorBehavior::Value::dummy);

    mod.add_type<autopas::options::IteratorBehavior>("IteratorBehavior")
            .constructor<autopas::options::IteratorBehavior::Value>();
    
    // add type MoleculeLJ
    mod.add_type<Parametric<TypeVar<1>>>("MoleculeLJ")
            .apply<autopas::MoleculeLJ<float>, autopas::MoleculeLJ<double>>(WrapMoleculeLJ());

    mod.add_type<Parametric<TypeVar<1>>>("ParticleBase")
            .apply<autopas::ParticleBase<float, int>, autopas::ParticleBase<double, long>>(WrapParticleBase());

    mod.add_type<Hydrogen>("Hydrogen")
            .constructor<>()
            .constructor<int, std::vector<double>, std::vector<double>>()
            .method("update_x", &Hydrogen::update_x)
            .method("update_v", &Hydrogen::update_v);

    mod.add_type<Parametric<TypeVar<1>>>("MoleculeJ")
            .apply<MoleculeJ<double>, MoleculeJ<float>>(WrapMoleculeJ()); 

    /*
    mod.add_type<autopas::options::IteratorBehavior>("IteratorBehavior", jlcxx::julia_base_type<autopas::options::Option<autopas::options::IteratorBehavior>>())
            .constructor<>();

    mod.add_type<Parametric<TypeVar<1>>>("Option")
            .apply<autopas::options::Option<autopas::options::IteratorBehavior>>(WrapOptions());

    */
    // add iterator types necessary to iterate over the AutoPasContainer
    // add type ParticleIteratorInterface
    mod.add_type<Parametric<TypeVar<1>, TypeVar<2>>>("ParticleIteratorInterface")
            .apply<autopas::ParticleIteratorInterface<autopas::MoleculeLJ<double>, true>>(WrapIteratorInterface());

    
    /*
    // add type ParticleIteratorInterfaceImpl
    mod.add_type<Parametric<TypeVar<1>, TypeVar<2>>>("ParticleIteratorInterfaceImpl", jlcxx::julia_base_type<autopas::ParticleIteratorInterface<autopas::MoleculeLJ<double>, true>>())
            .apply<autopas::internal::ParticleIteratorInterfaceImpl<autopas::MoleculeLJ<double>, true>>(WrapIteratorInterfaceImpl());
    */

    // add type ParticleIteratorWrapper which is used to iterate over the AutoPasContainer
    mod.add_type<Parametric<TypeVar<1>, TypeVar<2>>>("ParticleIteratorWrapper", jlcxx::julia_base_type<autopas::ParticleIteratorInterface<autopas::MoleculeLJ<double>, true>>())
            .apply<autopas::ParticleIteratorWrapper<autopas::MoleculeLJ<double>, true>>(WrapIteratorWrapper());
    
    // add type AutoPas
    // mod.add_type<Parametric<TypeVar<1>>>("AutoPas")
    //         .apply<autopas::AutoPas<autopas::ParticleBase<float, int>>, autopas::AutoPas<autopas::MoleculeLJ<double>>>(WrapAutoPas());
    
    mod.add_type<Parametric<TypeVar<1>>>("AutoPas")
            .apply<autopas::AutoPas<autopas::MoleculeLJ<double>>>(WrapAutoPas());

    /*
    mod.add_bits<ParticleType>("ParticleType", jlcxx::julia_type("CppEnum"));
    mod.set_const("owned", owned);
    */

    /*
     * adding all enums of the Optiontypes from AutoPas/options
     */
    
    /*
     * add enum IteratorBehavior
     */
    
    mod.add_bits<IteratorBehaviorJulia>("IteratorBehaviorBla", jlcxx::julia_type("CppEnum"));
    mod.set_const("owned1", owned1);
    mod.set_const("halo1", halo1);
    mod.set_const("ownedOrHalo1", ownedOrHalo1);
    mod.set_const("dummy1", dummy1);

    /*
     * add enum ContainerOptionJulia
     */
    mod.add_bits<autopas::options::ContainerOption::Value>("ContainerOptionValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("directSum", autopas::options::ContainerOption::Value::directSum);
    mod.set_const("linkedCells", autopas::options::ContainerOption::Value::linkedCells);
    mod.set_const("linkedCellsReferences", autopas::options::ContainerOption::Value::linkedCellsReferences);
    mod.set_const("varVerletListsAsBuild", autopas::options::ContainerOption::Value::varVerletListsAsBuild);
    mod.set_const("verletClusterLists", autopas::options::ContainerOption::Value::verletClusterLists);
    mod.set_const("verletLists", autopas::options::ContainerOption::Value::verletLists);
    mod.set_const("verletListsCells", autopas::options::ContainerOption::Value::verletListsCells);
    mod.set_const("pairwiseVerletLists", autopas::options::ContainerOption::Value::pairwiseVerletLists);
    mod.set_const("octree", autopas::options::ContainerOption::Value::octree);

    mod.add_type<autopas::options::ContainerOption>("ContainerOption")
            .constructor<autopas::options::ContainerOption::Value>();

    /*
     * add enum DataLayoutOptionJulia
     */
    mod.add_bits<autopas::options::DataLayoutOption::Value>("DataLayoutOptionValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("aos", autopas::options::DataLayoutOption::aos);
    mod.set_const("soa", autopas::options::DataLayoutOption::soa);

    /*
     * add enum Newton3OptionJulia
     */
    
    mod.add_bits<autopas::options::Newton3Option::Value>("Newton3OptionValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("disabled", autopas::options::Newton3Option::disabled);
    mod.set_const("enabled", autopas::options::Newton3Option::enabled);

    /*
     * add enum LoadEstimatorOptionJulia
     */
    mod.add_bits<autopas::LoadEstimatorOption::Value>("LoadEstimatorOptionValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("none", autopas::LoadEstimatorOption::none);
    mod.set_const("squaredParticlesPerCell", autopas::LoadEstimatorOption::squaredParticlesPerCell);
    mod.set_const("neighborListLength", autopas::LoadEstimatorOption::neighborListLength);

    /*
     * add enum TraversalOptionJulia
     */
    mod.add_bits<autopas::options::TraversalOption::Value>("TraversalOptionValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("ds_sequential", autopas::options::TraversalOption::ds_sequential);
    mod.set_const("lc_c01", autopas::options::TraversalOption::lc_c01);
    mod.set_const("lc_c01_combined_SoA", autopas::options::TraversalOption::lc_c01_combined_SoA);

    mod.method("setAllowedContainers", &setAllowedContainers);


    /*
     * add setters of AutoPas attributes which cannot directly be wrapped
     */
    
    /*
     * add setAllowedContainers
     */
    // mod.method("setAllowedContainers", &setAllowedContainers);

    /*
     * add setAllowedDataLayouts
     */
    // mod.method("setAllowedDataLayouts", setAllowedDataLayouts);

    /*
     * add setAllowedNewton3Options
     */
    // mod.method("setAllowedNewton3Options", setAllowedNewton3Options);

    /*
     * add setAllowedTraversals
     */
    // mod.method("setAllowedTraversals", setAllowedTraversals);

    /*
     * add setAllowedLoadEstimators
     */
    // mod.method("setAllowedLoadEstimators", setAllowedLoadEstimators);

    /*
     * add bein function of AutoPas
     */
    mod.method("begin", &begin_new);

    // mod.add_type<Parametric<TypeVar<1>>>("AutoPas")
    //         .apply<autopas::AutoPas<Hydrogen>, autopas::AutoPas<ParticleType>>(WrapAutoPas());

    mod.method("get_vec", &get_vec);
    mod.method("setBox", &setBox);
    mod.method("getBox", &getBox);
    mod.method("getContainerOption", &getContainerOption);
}

// cmake -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DAUTOPAS_BUILD_TESTS=OFF ..