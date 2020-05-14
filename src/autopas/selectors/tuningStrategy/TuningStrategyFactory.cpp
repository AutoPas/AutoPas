/**
 * @file TuningStrategyFactory.cpp
 * @author seckler
 * @date 07.02.20
 */

#include "TuningStrategyFactory.h"

#include "ActiveHarmony.h"
#include "BayesianSearch.h"
#include "FullSearch.h"
#include "FullSearchMPI.h"
#include "RandomSearch.h"

std::unique_ptr<autopas::TuningStrategyInterface> autopas::TuningStrategyFactory::generateTuningStrategy(
    autopas::TuningStrategyOption tuningStrategyOption, const std::set<autopas::ContainerOption> &allowedContainers,
    autopas::NumberSet<double> &allowedCellSizeFactors, const std::set<autopas::TraversalOption> &allowedTraversals,
    const std::set<autopas::DataLayoutOption> &allowedDataLayouts,
    const std::set<autopas::Newton3Option> &allowedNewton3Options, unsigned int maxEvidence,
    AcquisitionFunctionOption acquisitionFunctionOption, MPI_Comm comm) {
  // clang compiler bug requires static cast
  switch (static_cast<TuningStrategyOption>(tuningStrategyOption)) {
    case TuningStrategyOption::randomSearch: {
      return std::make_unique<RandomSearch>(allowedContainers, allowedCellSizeFactors, allowedTraversals,
                                            allowedDataLayouts, allowedNewton3Options, maxEvidence);
    }
    case TuningStrategyOption::fullSearch: {
      if (not allowedCellSizeFactors.isFinite()) {
        autopas::utils::ExceptionHandler::exception(
            "AutoPas::generateTuningStrategy: fullSearch can not handle infinite cellSizeFactors!");
        return nullptr;
      }

      return std::make_unique<FullSearch>(allowedContainers, allowedCellSizeFactors.getAll(), allowedTraversals,
                                          allowedDataLayouts, allowedNewton3Options);
    }

    case TuningStrategyOption::fullSearchMPI: {
      if (not allowedCellSizeFactors.isFinite()) {
        autopas::utils::ExceptionHandler::exception(
            "AutoPass::generateTuningStrategy: fullSearchMPI cannot handle infinite cellSizeFactors!");
        return nullptr;
      }

      return std::make_unique<FullSearchMPI>(allowedContainers, allowedCellSizeFactors.getAll(), allowedTraversals,
                                             allowedDataLayouts, allowedNewton3Options, comm);
    }

    case TuningStrategyOption::bayesianSearch: {
      return std::make_unique<BayesianSearch>(allowedContainers, allowedCellSizeFactors, allowedTraversals,
                                              allowedDataLayouts, allowedNewton3Options, maxEvidence,
                                              acquisitionFunctionOption);
    }

    case TuningStrategyOption::activeHarmony: {
      return std::make_unique<ActiveHarmony>(allowedContainers, allowedCellSizeFactors, allowedTraversals,
                                             allowedDataLayouts, allowedNewton3Options);
    }
  }

  autopas::utils::ExceptionHandler::exception("AutoPas::generateTuningStrategy: Unknown tuning strategy {}!",
                                              tuningStrategyOption);
  return nullptr;
}
