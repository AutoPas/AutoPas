/**
 * @file PredictionLogger.h
 * @author F. Gratl
 * @date 25/01/2021
 */

#pragma once

#include <spdlog/async.h>

#include "autopas/selectors/Configuration.h"

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
   */
  PredictionLogger();

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
                         const std::unordered_map<Configuration, size_t, ConfigHash> &configurationPredictions,
                         size_t predictionErrorValue, size_t tuningPhaseCounter);

 private:
  static inline std::string loggerName() { return "PredictionLogger"; };
};

}  // namespace autopas