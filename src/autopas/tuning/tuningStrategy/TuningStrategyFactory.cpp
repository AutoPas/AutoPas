/**
 * @file TuningStrategyFactory.cpp
 * @author seckler
 * @date 07.02.2020
 */

#include "TuningStrategyFactory.h"

#include "autopas/options/TuningStrategyOption.h"
#include "autopas/tuning/tuningStrategy/ActiveHarmony.h"
#include "autopas/tuning/tuningStrategy/BayesianClusterSearch.h"
#include "autopas/tuning/tuningStrategy/BayesianSearch.h"
#include "autopas/tuning/tuningStrategy/DeepReinforcementLearning.h"
#include "autopas/tuning/tuningStrategy/MPIParallelizedStrategy.h"
#include "autopas/tuning/tuningStrategy/PredictiveTuning.h"
#include "autopas/tuning/tuningStrategy/RandomSearch.h"
#include "autopas/tuning/tuningStrategy/ReinforcementLearning.h"
#include "autopas/tuning/tuningStrategy/SlowConfigFilter.h"
#include "autopas/tuning/tuningStrategy/SortByName.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyFactoryInfo.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/FuzzyTuning.h"
#include "autopas/tuning/tuningStrategy/ruleBasedTuning/RuleBasedTuning.h"
#include "autopas/tuning/utils/SearchSpaceGenerators.h"
#include "autopas/utils/NumberSetFinite.h"

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
                                                                const TuningStrategyOption tuningStrategyOption,
                                                                const TuningStrategyFactoryInfo &info,
                                                                const InteractionTypeOption interactionType,
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
          interactionType, searchSpaceDimensions.containerOptions,
          NumberSetFinite<double>{searchSpaceDimensions.cellSizeFactors}, searchSpaceDimensions.traversalOptions,
          searchSpaceDimensions.loadEstimatorOptions, searchSpaceDimensions.dataLayoutOptions,
          searchSpaceDimensions.newton3Options, info.maxEvidence, info.acquisitionFunctionOption);
      break;
    }

    case TuningStrategyOption::bayesianClusterSearch: {
      const auto searchSpaceDimensions = inferOptionDimensions(searchSpace);
      tuningStrategy = std::make_unique<BayesianClusterSearch>(
          interactionType, searchSpaceDimensions.containerOptions,
          NumberSetFinite<double>{searchSpaceDimensions.cellSizeFactors}, searchSpaceDimensions.traversalOptions,
          searchSpaceDimensions.loadEstimatorOptions, searchSpaceDimensions.dataLayoutOptions,
          searchSpaceDimensions.newton3Options, info.maxEvidence, info.acquisitionFunctionOption, outputSuffix);
      break;
    }

    case TuningStrategyOption::activeHarmony: {
      // If a AH-server is provided, but MPI is disallowed, we have to ignore the server.
      if (std::getenv("HARMONY_HOST") != nullptr and not info.mpiDivideAndConquer) {
        unsetenv("HARMONY_HOST");
        AutoPasLog(WARN,
                   "HARMONY_HOST is set to a value, but the MPI strategy option is set to noMPI. "
                   "HARMONY_HOST will be unset to enforce a local tuning session");
      }
      const auto searchSpaceDimensions = inferOptionDimensions(searchSpace);
      tuningStrategy = std::make_unique<ActiveHarmony>(
          interactionType, searchSpaceDimensions.containerOptions,
          NumberSetFinite<double>{searchSpaceDimensions.cellSizeFactors}, searchSpaceDimensions.traversalOptions,
          searchSpaceDimensions.loadEstimatorOptions, searchSpaceDimensions.dataLayoutOptions,
          searchSpaceDimensions.newton3Options, info.mpiDivideAndConquer, info.autopasMpiCommunicator);
      break;
    }

    case TuningStrategyOption::predictiveTuning: {
      tuningStrategy =
          std::make_unique<PredictiveTuning>(info.relativeOptimumRange, info.maxTuningPhasesWithoutTest,
                                             info.minNumberOfEvidence, info.extrapolationMethodOption, outputSuffix);
      break;
    }

    case TuningStrategyOption::reinforcementLearning: {
      tuningStrategy = std::make_unique<ReinforcementLearning>(searchSpace, info.learningRate, info.discountFactor);
      break;
    }

    case TuningStrategyOption::deepReinforcementLearning: {
      tuningStrategy =
          std::make_unique<DeepReinforcementLearning>(info.doReinforcementUpdates, info.numExplorationSamples);
      break;
    }

    case TuningStrategyOption::ruleBasedTuning: {
      tuningStrategy = std::make_unique<RuleBasedTuning>(searchSpace, /*verify mode*/ false, info.ruleFileName);
      break;
    }

    case TuningStrategyOption::fuzzyTuning: {
      tuningStrategy = std::make_unique<FuzzyTuning>(info.fuzzyRuleFileName);
      break;
    }

    case TuningStrategyOption::mpiDivideAndConquer: {
#ifndef AUTOPAS_INTERNODE_TUNING
      utils::ExceptionHandler::exception(
          "Cannot use the TuningStrategy mpiDivideAndConquer without AUTOPAS_INTERNODE_TUNING=ON.");
#endif
      if (not info.mpiDivideAndConquer) {
        AutoPasLog(WARN,
                   "Using TuningStrategy::mpiDivideAndConquer, but TuningStrategyFactoryInfo.mpiDivideAndConquer "
                   "is false. Other strategies will not be notified that the search space is split which might cause "
                   "problems.");
      }
      tuningStrategy = std::make_unique<MPIParallelizedStrategy>(
          MPIParallelizedStrategy::createFallBackConfiguration(searchSpace, interactionType),
          info.autopasMpiCommunicator, info.mpiTuningMaxDifferenceForBucket, info.mpiTuningWeightForMaxDensity);
      break;
    }

    case TuningStrategyOption::slowConfigFilter: {
      tuningStrategy = std::make_unique<SlowConfigFilter>(info.relativeBlacklistRange);
      break;
    }

    case TuningStrategyOption::sortByName: {
      tuningStrategy = std::make_unique<SortByName>();
      break;
    }

    default: {
      utils::ExceptionHandler::exception("AutoPas::generateTuningStrategy: Unknown tuning strategy {}!",
                                         tuningStrategyOption);
      break;
    }
  }
  return tuningStrategy;
}
}  // namespace autopas::TuningStrategyFactory