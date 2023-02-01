#pragma once
#include <jlcxx/jlcxx.hpp>

/**
 * define Julia module wtih AutoPas option types and enum values
 */
JLCXX_MODULE define_module_options(jlcxx::Module& mod) {

    /**
     * add AcquisitionFunctionOption values to Julia
     */
    mod.add_bits<autopas::options::AcquisitionFunctionOption::Value>("AcquisitionFunctionOptionValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("upperConfidenceBound", autopas::options::AcquisitionFunctionOption::Value::upperConfidenceBound);
    mod.set_const("mean", autopas::options::AcquisitionFunctionOption::Value::mean);
    mod.set_const("variance", autopas::options::AcquisitionFunctionOption::Value::variance);
    mod.set_const("probabilityOfImprovement", autopas::options::AcquisitionFunctionOption::Value::probabilityOfImprovement);
    mod.set_const("expectedImprovement", autopas::options::AcquisitionFunctionOption::Value::expectedImprovement);

    /**
     * add AcquisitionFunctionOption type and constructor to Julia
     */
    mod.add_type<autopas::options::AcquisitionFunctionOption>("AcquisitionFunctionOption")
            .constructor<autopas::options::AcquisitionFunctionOption::Value>();

    /**
     * add MPIStrategyOption values to Julia
     */
    mod.add_bits<autopas::options::MPIStrategyOption::Value>("MPIStrategyOptionValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("noMPI", autopas::options::MPIStrategyOption::Value::noMPI);
    mod.set_const("divideAndConquer", autopas::options::MPIStrategyOption::Value::divideAndConquer);

    /**
     * add MPIStrategyOption type and constructor to Julia
     */
    mod.add_type<autopas::options::MPIStrategyOption>("MPIStrategyOption")
            .constructor<autopas::options::MPIStrategyOption::Value>();

    /**
     * add TuningStrategyOption values to Julia
     */
    mod.add_bits<autopas::options::TuningStrategyOption::Value>("TuningStrategyOptionValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("randomSearch", autopas::options::TuningStrategyOption::Value::randomSearch);
    mod.set_const("fullSearch", autopas::options::TuningStrategyOption::Value::fullSearch);
    mod.set_const("bayesianSearch", autopas::options::TuningStrategyOption::Value::bayesianSearch);
    mod.set_const("bayesianClusterSearch", autopas::options::TuningStrategyOption::Value::bayesianClusterSearch);
    mod.set_const("activeHarmony", autopas::options::TuningStrategyOption::Value::activeHarmony);
    mod.set_const("predictiveTuning", autopas::options::TuningStrategyOption::Value::predictiveTuning);
    
    /**
     * add TuningStrategyOption type and constructor to Julia
     */
    mod.add_type<autopas::options::TuningStrategyOption>("TuningStrategyOption")
            .constructor<autopas::options::TuningStrategyOption::Value>();

    /**
     * add SelectorStrategyOption values to Julia
     */
    mod.add_bits<autopas::options::SelectorStrategyOption::Value>("SelectorStrategyOptionValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("fastestAbs", autopas::options::SelectorStrategyOption::Value::fastestAbs);
    mod.set_const("fastestMean", autopas::options::SelectorStrategyOption::Value::fastestMean);
    mod.set_const("fastestMedian", autopas::options::SelectorStrategyOption::Value::fastestMedian);

    /**
     * add SelectorStrategyOption type and constructor to Julia
     */
    mod.add_type<autopas::options::SelectorStrategyOption>("SelectorStrategyOption")
            .constructor<autopas::options::SelectorStrategyOption::Value>();

    /**
     * add ExtrapolationMethodOption values to Julia
     */
    mod.add_bits<autopas::options::ExtrapolationMethodOption::Value>("ExtrapolationMethodOptionValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("linePrediction", autopas::options::ExtrapolationMethodOption::Value::linePrediction);
    mod.set_const("linearRegression", autopas::options::ExtrapolationMethodOption::Value::linearRegression);
    mod.set_const("newton", autopas::options::ExtrapolationMethodOption::Value::newton);

    /**
     * add ExtrapolationMethodOption type and constructor to Julia
     */
    mod.add_type<autopas::options::ExtrapolationMethodOption>("ExtrapolationMethodOption")
            .constructor<autopas::options::ExtrapolationMethodOption::Value>();

    /**
     * add IteratorBehavior values to Julia
     */
    mod.add_bits<autopas::options::IteratorBehavior::Value>("IteratorBehaviorValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("owned", autopas::options::IteratorBehavior::Value::owned);
    mod.set_const("halo", autopas::options::IteratorBehavior::Value::halo);
    mod.set_const("ownedOrHalo", autopas::options::IteratorBehavior::Value::ownedOrHalo);
    mod.set_const("dummy", autopas::options::IteratorBehavior::Value::dummy);
    mod.set_const("ownedOrHaloOrDummy", autopas::options::IteratorBehavior::Value::ownedOrHaloOrDummy);
    mod.set_const("forceSequential", autopas::options::IteratorBehavior::Value::forceSequential);

    /**
     * add IteratorBehavior type and constructor to Julia
     */
    mod.add_type<autopas::options::IteratorBehavior>("IteratorBehavior")
            .constructor<autopas::options::IteratorBehavior::Value>();

    /*
     * add ContainerOption values to Julia
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

    /**
     * add ContainerOption type and constructor to Julia
     */
    mod.add_type<autopas::options::ContainerOption>("ContainerOption")
            .constructor<autopas::options::ContainerOption::Value>();

    /*
     * add DataLayoutOption values to Julia
     */
    mod.add_bits<autopas::options::DataLayoutOption::Value>("DataLayoutOptionValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("aos", autopas::options::DataLayoutOption::aos);
    mod.set_const("soa", autopas::options::DataLayoutOption::soa);

    /**
     * add DataLayoutOption type and constructor to Julia
     */
    mod.add_type<autopas::options::DataLayoutOption>("DataLayoutOption")
            .constructor<autopas::options::DataLayoutOption::Value>();

    /*
     * add Newton3Option values to Julia
     */
    mod.add_bits<autopas::options::Newton3Option::Value>("Newton3OptionValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("disabled", autopas::options::Newton3Option::disabled);
    mod.set_const("enabled", autopas::options::Newton3Option::enabled);

    /**
     * add Newton3Option type and constructor to Julia
     */
    mod.add_type<autopas::options::Newton3Option>("Newton3Option")
            .constructor<autopas::options::Newton3Option::Value>();

    /*
     * add LoadEstimatorOption values to Julia
     */
    mod.add_bits<autopas::LoadEstimatorOption::Value>("LoadEstimatorOptionValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("none", autopas::LoadEstimatorOption::none);
    mod.set_const("squaredParticlesPerCell", autopas::LoadEstimatorOption::squaredParticlesPerCell);
    mod.set_const("neighborListLength", autopas::LoadEstimatorOption::neighborListLength);

    /**
     * add LoadEstimatorOption type and constructor to Julia
     */
    mod.add_type<autopas::LoadEstimatorOption>("LoadEstimatorOption")
            .constructor<autopas::LoadEstimatorOption::Value>();

    /*
     * add TraversalOption values to Julia
     */
    mod.add_bits<autopas::options::TraversalOption::Value>("TraversalOptionValue", jlcxx::julia_type("CppEnum"));
    mod.set_const("ds_sequential", autopas::options::TraversalOption::ds_sequential);
    mod.set_const("lc_c01", autopas::options::TraversalOption::lc_c01);
    mod.set_const("lc_c01_combined_SoA", autopas::options::TraversalOption::lc_c01_combined_SoA);

    /**
     * add TraversalOption type and constructor to Julia
     */
    mod.add_type<autopas::options::TraversalOption>("TraversalOption")
            .constructor<autopas::options::TraversalOption::Value>();
}