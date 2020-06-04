/**
 * @file TuningStrategyFactory.cpp
 * @author seckler
 * @date 07.02.2020
 */

#include "TuningStrategyFactory.h"

#include "ActiveHarmony.h"
#include "BayesianClusterSearch.h"
#include "BayesianSearch.h"
#include "FullSearch.h"
#include "PredictiveTuning.h"
#include "RandomSearch.h"

std::unique_ptr<autopas::TuningStrategyInterface> autopas::TuningStrategyFactory::generateTuningStrategy(
    autopas::TuningStrategyOption tuningStrategyOption, const std::set<autopas::ContainerOption> &allowedContainers,
    autopas::NumberSet<double> &allowedCellSizeFactors, const std::set<autopas::TraversalOption> &allowedTraversals,
    const std::set<autopas::LoadEstimatorOption> &allowedLoadEstimators,
    const std::set<autopas::DataLayoutOption> &allowedDataLayouts,
    const std::set<autopas::Newton3Option> &allowedNewton3Options, unsigned int maxEvidence, double relativeOptimum,
    unsigned int maxTuningPhasesWithoutTest, AcquisitionFunctionOption acquisitionFunctionOption) {
  // clang compiler bug requires static cast
  switch (static_cast<TuningStrategyOption>(tuningStrategyOption)) {
    case TuningStrategyOption::randomSearch: {
      return std::make_unique<RandomSearch>(allowedContainers, allowedCellSizeFactors, allowedTraversals,
                                            allowedLoadEstimators, allowedDataLayouts, allowedNewton3Options,
                                            maxEvidence);
    }
    case TuningStrategyOption::fullSearch: {
      if (not allowedCellSizeFactors.isFinite()) {
        autopas::utils::ExceptionHandler::exception(
            "AutoPas::generateTuningStrategy: fullSearch can not handle infinite cellSizeFactors!");
        return nullptr;
      }

      return std::make_unique<FullSearch>(allowedContainers, allowedCellSizeFactors.getAll(), allowedTraversals,
                                          allowedLoadEstimators, allowedDataLayouts, allowedNewton3Options);
    }

    case TuningStrategyOption::bayesianSearch: {
      return std::make_unique<BayesianSearch>(allowedContainers, allowedCellSizeFactors, allowedTraversals,
                                              allowedLoadEstimators, allowedDataLayouts, allowedNewton3Options,
                                              maxEvidence, acquisitionFunctionOption);
    }

    case TuningStrategyOption::bayesianClusterSearch: {
      return std::make_unique<BayesianClusterSearch>(allowedContainers, allowedCellSizeFactors, allowedTraversals,
                                                     allowedLoadEstimators, allowedDataLayouts, allowedNewton3Options,
                                                     maxEvidence, acquisitionFunctionOption);
    }

    case TuningStrategyOption::activeHarmony: {
      return std::make_unique<ActiveHarmony>(allowedContainers, allowedCellSizeFactors, allowedTraversals,
                                             allowedLoadEstimators, allowedDataLayouts, allowedNewton3Options);
    }

    case TuningStrategyOption::predictiveTuning: {
      return std::make_unique<PredictiveTuning>(allowedContainers, allowedCellSizeFactors.getAll(), allowedTraversals,
                                                allowedLoadEstimators, allowedDataLayouts, allowedNewton3Options,
                                                relativeOptimum, maxTuningPhasesWithoutTest);
    }
  }

  autopas::utils::ExceptionHandler::exception("AutoPas::generateTuningStrategy: Unknown tuning strategy {}!",
                                              tuningStrategyOption);
  return nullptr;
}
