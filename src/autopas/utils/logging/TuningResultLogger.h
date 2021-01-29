/**
 * @file TuningResultLogger.h
 * @author F. Gratl
 * @date 25/01/2021
 */

#pragma once

#include <spdlog/async.h>

#include "autopas/selectors/Configuration.h"

namespace autopas {

/**
 * Helper to log results of the tuning process to a csv file for easier analysis.
 *
 * It uses an asynchronous spd logger to write a csv file named "AutoPas_iterationPerformance_<dateStamp>.csv".
 *
 * By default logging the data is disabled. It can be enabled by setting the cmake variable AUTOPAS_LOG_TUNINGRESULTS
 * to ON.
 */
class TuningResultLogger {
 public:
  /**
   * Constructor initializes the logger and sets the output file name.
   */
  explicit TuningResultLogger(const std::string &outputSuffix = "");

  /**
   * Destructor drops the logger from the spd registry.
   */
  ~TuningResultLogger();

  /**
   * Log the result of a tuning phase.
   * @param configuration
   * @param iteration
   * @param timeTuning
   */
  void logTuningResult(const autopas::Configuration &configuration, size_t iteration, long timeTuning);

 private:
  static inline std::string loggerName() { return "TuningResultLogger"; };
};

}  // namespace autopas