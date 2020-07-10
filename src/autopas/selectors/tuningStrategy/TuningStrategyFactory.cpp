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
#include "autopas/options/MPIStrategyOption.h"
#include "autopas/utils/AutoPasConfigurationCommunicator.h"

std::unique_ptr<autopas::TuningStrategyInterface> autopas::TuningStrategyFactory::generateTuningStrategy(
    autopas::TuningStrategyOption tuningStrategyOption, std::set<autopas::ContainerOption> &allowedContainers,
    autopas::NumberSet<double> &allowedCellSizeFactors, std::set<autopas::TraversalOption> &allowedTraversals,
    std::set<autopas::LoadEstimatorOption> &allowedLoadEstimators,
    std::set<autopas::DataLayoutOption> &allowedDataLayouts, std::set<autopas::Newton3Option> &allowedNewton3Options,
    unsigned int maxEvidence, double relativeOptimum, unsigned int maxTuningPhasesWithoutTest,
    unsigned int evidenceFirstPrediction, AcquisitionFunctionOption acquisitionFunctionOption,
    ExtrapolationMethodOption extrapolationMethodOption, MPIStrategyOption mpiStrategyOption, AutoPas_MPI_Comm comm) {
  switch (static_cast<autopas::MPIStrategyOption>(mpiStrategyOption)) {
    case MPIStrategyOption::noMPI: {
      break;
    }

    case MPIStrategyOption::divideAndConquer: {
      if (tuningStrategyOption == TuningStrategyOption::activeHarmony && getenv("HARMONY_HOST") != nullptr) {
        // rank 0 will solely set up the entire search, so we cannot divide the search space
        break;
      }
      int rank, commSize;
      AutoPas_MPI_Comm_rank(comm, &rank);
      AutoPas_MPI_Comm_size(comm, &commSize);
      utils::AutoPasConfigurationCommunicator::distributeConfigurations(
          allowedContainers, allowedCellSizeFactors, allowedTraversals, allowedLoadEstimators, allowedDataLayouts,
          allowedNewton3Options, rank, commSize);
      break;
    }

    default: {
      utils::ExceptionHandler::exception("AutoPas::generateTuningStrategy: Unknown MPIStrategyOption: {}",
                                         mpiStrategyOption);
    }
  }

  std::unique_ptr<autopas::TuningStrategyInterface> tuningStrategy = nullptr;
  switch (static_cast<autopas::TuningStrategyOption>(tuningStrategyOption)) {
    case TuningStrategyOption::randomSearch: {
      tuningStrategy =
          std::make_unique<RandomSearch>(allowedContainers, allowedCellSizeFactors, allowedTraversals,
                                         allowedLoadEstimators, allowedDataLayouts, allowedNewton3Options, maxEvidence);
      break;
    }

    case TuningStrategyOption::fullSearch: {
      if (not allowedCellSizeFactors.isFinite()) {
        autopas::utils::ExceptionHandler::exception(
            "AutoPas::generateTuningStrategy: fullSearch can not handle infinite cellSizeFactors!");
        return nullptr;
      }

      tuningStrategy =
          std::make_unique<FullSearch>(allowedContainers, allowedCellSizeFactors.getAll(), allowedTraversals,
                                       allowedLoadEstimators, allowedDataLayouts, allowedNewton3Options);
      break;
    }

    case TuningStrategyOption::bayesianSearch: {
      tuningStrategy = std::make_unique<BayesianSearch>(allowedContainers, allowedCellSizeFactors, allowedTraversals,
                                                        allowedLoadEstimators, allowedDataLayouts,
                                                        allowedNewton3Options, maxEvidence, acquisitionFunctionOption);
      break;
    }

    case TuningStrategyOption::bayesianClusterSearch: {
      tuningStrategy = std::make_unique<BayesianClusterSearch>(
          allowedContainers, allowedCellSizeFactors, allowedTraversals, allowedLoadEstimators, allowedDataLayouts,
          allowedNewton3Options, maxEvidence, acquisitionFunctionOption);
      break;
    }

    case TuningStrategyOption::activeHarmony: {
      // If a AH-server is provided, but MPI is disallowed, we have to ignore the server
      if (getenv("HARMONY_HOST") != nullptr && mpiStrategyOption == MPIStrategyOption::noMPI) {
        unsetenv("HARMONY_HOST");
        AutoPasLog(warn,
                   "HARMONY_HOST is set to a value, but the MPI strategy option is set to noMPI. "
                   "HARMONY_HOST will be unset to enforce a local tuning session");
      }
      tuningStrategy = std::make_unique<ActiveHarmony>(allowedContainers, allowedCellSizeFactors, allowedTraversals,
                                                       allowedLoadEstimators, allowedDataLayouts, allowedNewton3Options,
                                                       mpiStrategyOption, comm);
      break;
    }

    case TuningStrategyOption::predictiveTuning: {
      tuningStrategy = std::make_unique<PredictiveTuning>(
          allowedContainers, allowedCellSizeFactors.getAll(), allowedTraversals, allowedLoadEstimators,
          allowedDataLayouts, allowedNewton3Options, relativeOptimum, maxTuningPhasesWithoutTest,
          evidenceFirstPrediction, extrapolationMethodOption);
      break;
    }

    default: {
      autopas::utils::ExceptionHandler::exception("AutoPas::generateTuningStrategy: Unknown tuning strategy {}!",
                                                  tuningStrategyOption);
      break;
    }
  }

  switch (static_cast<MPIStrategyOption>(mpiStrategyOption)) {
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
      return std::make_unique<MPIParallelizedStrategy>(std::move(tuningStrategy), comm);
    }
  }

  return nullptr;
}
