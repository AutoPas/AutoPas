/**
 * @file TuningDataLogger.h
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
 * By default logging the data is disabled. It can be enabled by setting the cmake variable AUTOPAS_LOG_TUNINGDATA
 * to ON.
 */
class TuningDataLogger {
 public:
  /**
   * Constructor initializes the logger and sets the output file name.
   * @param numSamples Number of samples that are taken per configuration.
   * @param outputSuffix Suffix for all output files produced by this class.
   */
  explicit TuningDataLogger(size_t numSamples, const std::string &outputSuffix = "");

  /**
   * Destructor drops the logger from the spd registry.
   */
  ~TuningDataLogger();

  /**
   * Log the result of a tuning phase.
   * @param configuration
   * @param samples
   * @param iteration
   * @param reducedValue
   * @param smoothedVale
   */
  void logTuningData(const autopas::Configuration &configuration, const std::vector<size_t> &samples, size_t iteration,
                     size_t reducedValue, size_t smoothedVale);

 private:
  static inline std::string loggerName() { return "TuningDataLogger"; };
};

}  // namespace autopas