/**
 * @file PredictionLogger.h
 * @author F. Gratl
 * @date 25/01/2021
 */

#pragma once

#include <spdlog/async.h>

#include "autopas/tuning/Configuration.h"

namespace autopas {

/**
 * Helper to log prediction data of PredictiveTuning to a csv file for easier analysis.
 *
 * It uses an asynchronous spd logger to write a csv file named "AutoPas_predictions_<dateStamp>.csv".
 *
 * By default logging the data is disabled. It can be enabled by setting the cmake variable AUTOPAS_LOG_PREDICTIONS
 * to ON.
 */
class PredictionLogger {
 public:
  /**
   * Constructor initializes the logger and sets the output file name.
   * @param outputSuffix Suffix for all output files produced by this class.
   */
  explicit PredictionLogger(const std::string &outputSuffix = "");

  /**
   * Destructor drops the logger from the spd registry.
   */
  ~PredictionLogger();

  /**
   * Print all predictions of the given set of configurations to the logger.
   * @param configurations
   * @param configurationPredictions
   * @param predictionErrorValue
   * @param tuningPhaseCounter
   */
  void logAllPredictions(const std::set<Configuration> &configurations,
                         const std::unordered_map<Configuration, long, ConfigHash> &configurationPredictions,
                         long predictionErrorValue, size_t tuningPhaseCounter);

 private:
  std::string _loggerName;
};

}  // namespace autopas