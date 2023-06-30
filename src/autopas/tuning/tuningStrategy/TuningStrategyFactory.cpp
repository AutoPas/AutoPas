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
#include "MPIParallelizedStrategy.h"
#include "PredictiveTuning.h"
#include "RandomSearch.h"
#include "TuningStrategyFactoryInfo.h"
#include "autopas/options/MPIStrategyOption.h"
#include "autopas/tuning/tuningStrategy/ruleBasedTuning/RuleBasedTuning.h"
#include "autopas/utils/AutoPasConfigurationCommunicator.h"

std::unique_ptr<autopas::TuningStrategyInterface> autopas::TuningStrategyFactory::generateTuningStrategy(
    const std::set<Configuration> &searchSpace, autopas::TuningStrategyOption tuningStrategyOption,
    const TuningStrategyFactoryInfo &info, const std::string &outputSuffix) {
  // ======== prepare MPI =====================================================

  // only needed in the MPI case, but need to be declared here.
  std::set<autopas::ContainerOption> fallbackContainers;
  std::unique_ptr<autopas::NumberSet<double>> fallbackCellSizeFactors;
  std::set<autopas::TraversalOption> fallbackTraversals;
  std::set<autopas::LoadEstimatorOption> fallbackLoadEstimators;
  std::set<autopas::DataLayoutOption> fallbackDataLayouts;
  std::set<autopas::Newton3Option> fallbackNewton3;
  // if an mpi-strategy is used, the local search space is set up here, as well as the fallback options.
  switch (static_cast<autopas::MPIStrategyOption>(info.mpiStrategyOption)) {
    case MPIStrategyOption::noMPI: {
      break;
    }

    case MPIStrategyOption::divideAndConquer: {
#ifndef AUTOPAS_INTERNODE_TUNING
      utils::ExceptionHandler::exception(
          "Cannot use the divideAndConquer search-strategy without AUTOPAS_INTERNODE_TUNING=ON."
          "aborting.");
#endif
      if (tuningStrategyOption == TuningStrategyOption::activeHarmony and getenv("HARMONY_HOST") != nullptr) {
        // rank 0 will solely set up the entire search, so we cannot divide the search space.
        break;
      }
      int rank, commSize;
      AutoPas_MPI_Comm_rank(info.comm, &rank);
      AutoPas_MPI_Comm_size(info.comm, &commSize);
      fallbackContainers = std::set<autopas::ContainerOption>(allowedContainers);
      if (allowedCellSizeFactors.isFinite()) {
        fallbackCellSizeFactors = std::make_unique<autopas::NumberSetFinite<double>>(allowedCellSizeFactors.getAll());
      } else {
        fallbackCellSizeFactors = std::make_unique<autopas::NumberInterval<double>>(allowedCellSizeFactors.getMin(),
                                                                                    allowedCellSizeFactors.getMax());
      }
      fallbackTraversals = std::set<autopas::TraversalOption>(allowedTraversals);
      fallbackLoadEstimators = std::set<autopas::LoadEstimatorOption>(allowedLoadEstimators);
      fallbackDataLayouts = std::set<autopas::DataLayoutOption>(allowedDataLayouts);
      fallbackNewton3 = std::set<autopas::Newton3Option>(allowedNewton3Options);

      utils::AutoPasConfigurationCommunicator::distributeConfigurations(
          allowedContainers, allowedCellSizeFactors, allowedTraversals, allowedLoadEstimators, allowedDataLayouts,
          allowedNewton3Options, rank, commSize);
      break;
    }

    default: {
      utils::ExceptionHandler::exception("AutoPas::generateTuningStrategy: Unknown MPIStrategyOption: {}",
                                         info.mpiStrategyOption);
    }
  }

  // ======== initiate tuning strategy ========================================

  std::unique_ptr<autopas::TuningStrategyInterface> tuningStrategy = nullptr;
  switch (static_cast<autopas::TuningStrategyOption>(tuningStrategyOption)) {
    case TuningStrategyOption::randomSearch: {
      tuningStrategy = std::make_unique<RandomSearch>(info.maxEvidence, std::random_device()());
      break;
    }

    case TuningStrategyOption::fullSearch: {
      autopas::utils::ExceptionHandler::exception(
          "Full Search is no tuning strategy anymore! If you want this behavior don't select any tuning strategy.");
      break;
    }

    case TuningStrategyOption::bayesianSearch: {
      tuningStrategy = std::make_unique<BayesianSearch>(
          allowedContainers, allowedCellSizeFactors, allowedTraversals, allowedLoadEstimators, allowedDataLayouts,
          allowedNewton3Options, info.maxEvidence, info.acquisitionFunctionOption);
      break;
    }

    case TuningStrategyOption::bayesianClusterSearch: {
      tuningStrategy = std::make_unique<BayesianClusterSearch>(
          allowedContainers, allowedCellSizeFactors, allowedTraversals, allowedLoadEstimators, allowedDataLayouts,
          allowedNewton3Options, info.maxEvidence, info.acquisitionFunctionOption, outputSuffix);
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
      tuningStrategy = std::make_unique<ActiveHarmony>(allowedContainers, allowedCellSizeFactors, allowedTraversals,
                                                       allowedLoadEstimators, allowedDataLayouts, allowedNewton3Options,
                                                       info.mpiStrategyOption, info.comm);
      break;
    }

    case TuningStrategyOption::predictiveTuning: {
      tuningStrategy =
          std::make_unique<PredictiveTuning>(searchSpace, info.relativeOptimum, info.maxTuningPhasesWithoutTest,
                                             info.minNumberOfEvidence, info.extrapolationMethodOption, outputSuffix);
      break;
    }

    case TuningStrategyOption::ruleBasedTuning: {
      if (not allowedCellSizeFactors.isFinite()) {
        // FIXME: This is probably not true anymore
        autopas::utils::ExceptionHandler::exception(
            "AutoPas::generateTuningStrategy: ruleBasedTuning can not handle infinite cellSizeFactors!");
        return nullptr;
      }

      tuningStrategy = std::make_unique<RuleBasedTuning>(
          allowedContainers, allowedCellSizeFactors.getAll(), allowedTraversals, allowedLoadEstimators,
          allowedDataLayouts, allowedNewton3Options, /*verify mode*/ false, info.ruleFileName);
      break;
    }

    default: {
      autopas::utils::ExceptionHandler::exception("AutoPas::generateTuningStrategy: Unknown tuning strategy {}!",
                                                  tuningStrategyOption);
      break;
    }
  }

  // ======== Wrap strategy into MPI wrapper if appropriate ===================

  switch (static_cast<MPIStrategyOption>(info.mpiStrategyOption)) {
    case MPIStrategyOption::noMPI: {
      return tuningStrategy;
    }

    case MPIStrategyOption::divideAndConquer: {
      if (tuningStrategyOption == TuningStrategyOption::activeHarmony) {
        if (getenv("HARMONY_HOST") != nullptr) {
          // A server has been specified, so no need to handle communication via MPI as well.
          return tuningStrategy;
        }
      }
      return std::make_unique<MPIParallelizedStrategy>(std::move(tuningStrategy), info.comm, fallbackContainers,
                                                       fallbackTraversals, fallbackLoadEstimators, fallbackDataLayouts,
                                                       fallbackNewton3);
    }
  }

  return nullptr;
}
