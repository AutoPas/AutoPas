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
    AcquisitionFunctionOption acquisitionFunctionOption, MPIStrategyOption mpiStrategyOption, AutoPas_MPI_Comm comm) {
  switch (static_cast<autopas::MPIStrategyOption>(mpiStrategyOption)) {
    case MPIStrategyOption::noMPI: {
      break;
    }

    case MPIStrategyOption::divideAndConquer: {
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
      tuningStrategy =
          std::make_unique<ActiveHarmony>(allowedContainers, allowedCellSizeFactors, allowedTraversals,
                                          allowedLoadEstimators, allowedDataLayouts, allowedNewton3Options);
      break;
    }

    case TuningStrategyOption::predictiveTuning: {
      tuningStrategy = std::make_unique<PredictiveTuning>(
          allowedContainers, allowedCellSizeFactors.getAll(), allowedTraversals, allowedLoadEstimators,
          allowedDataLayouts, allowedNewton3Options, relativeOptimum, maxTuningPhasesWithoutTest);
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
      return std::make_unique<MPIParallelizedStrategy>(std::move(tuningStrategy), comm);
    }
  }

  return nullptr;
}
