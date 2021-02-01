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
   * Log the tuning time to the loggers buffer.
   * @param timeTuning
   */
  void logTimeTuning(long timeTuning);

  /**
   * Log the given arguments and the internal buffer to the csv file.
   * @param configuration
   * @param iteration
   * @param inTuningPhase
   * @param timeIteratePairwise
   * @param timeRebuildNeighborLists
   * @param timeWholeIteration
   */
  void logIteration(const Configuration &configuration, size_t iteration, bool inTuningPhase, long timeIteratePairwise,
                    long timeRebuildNeighborLists, long timeWholeIteration);

 private:
  std::string _loggerName;
#ifdef AUTOPAS_LOG_ITERATIONS
  long _bufferTimeTuning{0};
#endif
};

}  // namespace autopas