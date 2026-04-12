/**
 * @file TuningResultLogger.h
 * @author F. Gratl
 * @date 25/01/2021
 */

#pragma once

#include <spdlog/async.h>

#include "autopas/options/TuningMetricOption.h"
#include "autopas/tuning/Configuration.h"

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
   * @param outputSuffix Suffix for all output files produced by this class.
   * @param tuningMetric Tuning metric (time or energy) used for the current simulation.
   */
  explicit TuningResultLogger(const std::string &outputSuffix = "",
                              TuningMetricOption tuningMetric = TuningMetricOption::time);

  /**
   * Destructor drops the logger from the spd registry.
   */
  ~TuningResultLogger();

  /**
   * Log the result of a tuning phase.
   * @param configuration
   * @param iteration
   * @param timeTuning
   * @param optimumPerformance Performance of the best configuration
   */
  void logTuningResult(const autopas::Configuration &configuration, size_t iteration, long timeTuning,
                       long optimumPerformance) const;

 private:
  std::string _loggerName;
};

}  // namespace autopas