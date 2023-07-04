/**
 * @file TuningStrategyFactory.cpp
 * @author seckler
 * @date 07.02.2020
 */

#include "TuningStrategyFactory.h"

#include "ActiveHarmony.h"
#include "BayesianClusterSearch.h"
#include "BayesianSearch.h"
#include "MPIParallelizedStrategy.h"
#include "PredictiveTuning.h"
#include "RandomSearch.h"
#include "TuningStrategyFactoryInfo.h"
#include "autopas/options/MPIStrategyOption.h"
#include "autopas/tuning/tuningStrategy/ruleBasedTuning/RuleBasedTuning.h"
#include "autopas/utils/AutoPasConfigurationCommunicator.h"
#include "options/TuningStrategyOption.h"
#include "tuning/utils/SearchSpaceGenerators.h"
#include "utils/NumberSetFinite.h"

namespace autopas::TuningStrategyFactory {

/**
 * Wraps SearchSpaceGenerators::inferOptionDimensions() and adds a warning about its usage.
 *
 * This function acts as a workaround for old and complex tuning strategies that rely on the search space to
 * be represented as a set of vectors of available options.
 *
 * @param searchSpace
 * @return
 */
SearchSpaceGenerators::OptionSpace inferOptionDimensions(const std::set<Configuration> &searchSpace) {
  AutoPasLog(WARN,
             "Inferring the dimensions of the search space from the given set of configurations."
             "This only works reliably if the set was created from a cross product of option vectors.");
  return SearchSpaceGenerators::inferOptionDimensions(searchSpace);
}

std::unique_ptr<TuningStrategyInterface> generateTuningStrategy(const std::set<Configuration> &searchSpace,
                                                                TuningStrategyOption tuningStrategyOption,
                                                                const TuningStrategyFactoryInfo &info,
                                                                const std::string &outputSuffix) {
  std::unique_ptr<TuningStrategyInterface> tuningStrategy = nullptr;
  switch (static_cast<TuningStrategyOption>(tuningStrategyOption)) {
    case TuningStrategyOption::randomSearch: {
      tuningStrategy = std::make_unique<RandomSearch>(info.maxEvidence, std::random_device()());
      break;
    }

    case TuningStrategyOption::fullSearch: {
      utils::ExceptionHandler::exception(
          "Full Search is no tuning strategy anymore! If you want this behavior don't select any tuning strategy.");
      break;
    }

    case TuningStrategyOption::bayesianSearch: {
      const auto searchSpaceDimensions = inferOptionDimensions(searchSpace);
      tuningStrategy = std::make_unique<BayesianSearch>(
          searchSpaceDimensions.containerOptions, NumberSetFinite<double>{searchSpaceDimensions.cellSizeFactors},
          searchSpaceDimensions.traversalOptions, searchSpaceDimensions.loadEstimatorOptions,
          searchSpaceDimensions.dataLayoutOptions, searchSpaceDimensions.newton3Options, info.maxEvidence,
          info.acquisitionFunctionOption);
      break;
    }

    case TuningStrategyOption::bayesianClusterSearch: {
      const auto searchSpaceDimensions = inferOptionDimensions(searchSpace);
      tuningStrategy = std::make_unique<BayesianClusterSearch>(
          searchSpaceDimensions.containerOptions, NumberSetFinite<double>{searchSpaceDimensions.cellSizeFactors},
          searchSpaceDimensions.traversalOptions, searchSpaceDimensions.loadEstimatorOptions,
          searchSpaceDimensions.dataLayoutOptions, searchSpaceDimensions.newton3Options, info.maxEvidence,
          info.acquisitionFunctionOption, outputSuffix);
      break;
    }

    case TuningStrategyOption::activeHarmony: {
      // If a AH-server is provided, but MPI is disallowed, we have to ignore the server.
      if (std::getenv("HARMONY_HOST") != nullptr and info.mpiStrategyOption == MPIStrategyOption::noMPI) {
        unsetenv("HARMONY_HOST");
        AutoPasLog(WARN,
                   "HARMONY_HOST is set to a value, but the MPI strategy option is set to noMPI. "
                   "HARMONY_HOST will be unset to enforce a local tuning session");
      }
      const auto searchSpaceDimensions = inferOptionDimensions(searchSpace);
      tuningStrategy = std::make_unique<ActiveHarmony>(
          searchSpaceDimensions.containerOptions, NumberSetFinite<double>{searchSpaceDimensions.cellSizeFactors},
          searchSpaceDimensions.traversalOptions, searchSpaceDimensions.loadEstimatorOptions,
          searchSpaceDimensions.dataLayoutOptions, searchSpaceDimensions.newton3Options, info.mpiStrategyOption,
          info.comm);
      break;
    }

    case TuningStrategyOption::predictiveTuning: {
      tuningStrategy =
          std::make_unique<PredictiveTuning>(info.relativeOptimum, info.maxTuningPhasesWithoutTest,
                                             info.minNumberOfEvidence, info.extrapolationMethodOption, outputSuffix);
      break;
    }

    case TuningStrategyOption::ruleBasedTuning: {
      tuningStrategy = std::make_unique<RuleBasedTuning>(searchSpace, /*verify mode*/ false, info.ruleFileName);
      break;
    }

    case TuningStrategyOption::mpiDivideAndConquer: {
#ifndef AUTOPAS_INTERNODE_TUNING
      utils::ExceptionHandler::exception(
          "Cannot use the TuningStrategy mpiDivideAndConquer without AUTOPAS_INTERNODE_TUNING=ON.");
#endif
      tuningStrategy = std::make_unique<MPIParallelizedStrategy>(
          MPIParallelizedStrategy::createFallBackConfiguration(searchSpace), info.comm,
          info.mpiTuningMaxDifferenceForBucket, info.mpiTuningWeightForMaxDensity);
      break;
    }

    default: {
      utils::ExceptionHandler::exception("AutoPas::generateTuningStrategy: Unknown tuning strategy {}!",
                                         tuningStrategyOption);
      break;
    }
  }
  return nullptr;
}
}  // namespace autopas::TuningStrategyFactory