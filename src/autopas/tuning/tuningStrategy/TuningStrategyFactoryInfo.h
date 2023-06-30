/**
 * @file TuningStrategyFactoryInfo.h
 * @author F. Gratl
 * @date 30.06.23
 */

#pragma once

#include <string>

#include "autopas/options/AcquisitionFunctionOption.h"
#include "autopas/options/ExtrapolationMethodOption.h"
#include "autopas/options/MPIStrategyOption.h"
#include "autopas/utils/WrapMPI.h"

namespace autopas {

/**
 * Helper struct encapsulating most information for the tuning factory.
 * This way if only a specific option should be built, not all arguments have to be specified explicitly.
 */
struct TuningStrategyFactoryInfo {
  // Used by multiple Strategies
  /**
   * If the strategy doesn't converge (or not quickly enough) limit the number of evidence to this number.
   */
  unsigned int maxEvidence{10};

  // Predictive Tuning Options
  /**
   * Function option used for extrapolating performance from observed evidence.
   */
  ExtrapolationMethodOption extrapolationMethodOption{ExtrapolationMethodOption::linearRegression};
  /**
   * Factor of the range of the optimal configurations for the optimalSearchSpace.
   */
  double relativeOptimum{1.2};
  /**
   * If a config is not tested for this number of tuning phases test it again to make predictions more reliable.
   */
  unsigned int maxTuningPhasesWithoutTest{100};
  /**
   * The number of evidence that have to be collected until the first prediction can be made.
   */
  unsigned int minNumberOfEvidence{3};

  // Slow Config Filter Options
  /**
   * Any configuration that is slower than the fastest times this factor will be blacklisted.
   */
  double relativeBlacklistRange{3};

  // Bayesian Strategies Options
  /**
   * Function used to predict informational gain.
   */
  AcquisitionFunctionOption acquisitionFunctionOption{AcquisitionFunctionOption::upperConfidenceBound};

  // Rule Based Tuning Options
  /**
   * The name and path of the file where the rules are stored.
   */
  std::string ruleFileName{};

  // MPI Tuning Options
  /**
   * Whether the chosen tuning strategy will be parallelized by the MPI wrapper.
   */
  MPIStrategyOption mpiStrategyOption = MPIStrategyOption::noMPI;
  /**
   * MPI Communicator used within AutoPas.
   */
  AutoPas_MPI_Comm comm = AUTOPAS_MPI_COMM_WORLD;
};
}  // namespace autopas
