/**
 * @file IterationLogger.h
 * @author F. Gratl
 * @date 25/01/2021
 */

#pragma once

#include <spdlog/async.h>

#include "autopas/selectors/Configuration.h"

namespace autopas {

/**
 * Helper to log performance data of AutoPas::iteratePairwise() to a csv file for easier analysis.
 *
 * It uses an asynchronous spd logger to write a csv file named "AutoPas_iterationPerformance_<dateStamp>.csv".
 *
 * By default logging the data is disabled. It can be enabled by setting the cmake variable AUTOPAS_LOG_ITERATIONS
 * to ON.
 */
class IterationLogger {
 public:
  /**
   * Constructor initializes the logger and sets the output file name.
   * @param outputSuffix Suffix for all output files produced by this class.
   */
  explicit IterationLogger(const std::string &outputSuffix = "");

  /**
   * Destructor drops the logger from the spd registry.
   */
  ~IterationLogger();

  /**
   * Log the given arguments and the internal buffer to the csv file.
   * @param configuration
   * @param iteration
   * @param inTuningPhase True if the logged iteration is in a tuning phase.
   * @param timeIteratePairwise Time for Container::iteratePairwise().
   * @param timeRemainderTraversal Time for LogicHandler::doRemainderTraversal().
   * @param timeRebuildNeighborLists Time for Container::rebuildNeighborLists().
   * @param timeIteratePairwiseTotal Time for LogicHandler::iteratePairwise(). This is slightly more than the sum of the
   * above. Additional steps, only included in this timer are e.g. Functor::initTraversal() and Functor::endTraversal().
   * @param timeTuning Time for finding the next configuration.
   * @param energyPsys Energy in Joules for the entire measured system since the last measurement.
   * @param energyPkg Energy in Joules for the CPU package since the last measurement.
   * @param energyRam Energy in Joules for the RAM since the last measurement.
   */
  void logIteration(const Configuration &configuration, size_t iteration, bool inTuningPhase, long timeIteratePairwise,
                    long timeRemainderTraversal, long timeRebuildNeighborLists, long timeIteratePairwiseTotal,
                    long timeTuning, double energyPsys, double energyPkg, double energyRam);

 private:
  std::string _loggerName;
};

}  // namespace autopas
