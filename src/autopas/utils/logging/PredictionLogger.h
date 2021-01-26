/**
 * @file PredictionLogger.h
 * @author F. Gratl
 * @date 25/01/2021
 */

#pragma once

#include <spdlog/async.h>

#include "selectors/Configuration.h"
#include "utils/Timer.h"
#include "utils/logging/Logger.h"

namespace autopas {

class PredictionLogger {
 public:
  PredictionLogger();

  ~PredictionLogger();

  /**
   * Print all predictions of the given set of configurations to the logger.
   * @param configurations
   * @param _configurationPredictions
   * @param _predictionErrorValue
   * @param _tuningPhaseCounter
   */
  void logAllPredictions(const std::set<Configuration> &configurations,
                         const std::unordered_map<Configuration, size_t, ConfigHash> &configurationPredictions,
                         size_t predictionErrorValue, size_t tuningPhaseCounter);

 private:
  static inline std::string loggerName() { return "PredictionLogger"; };
};

}  // namespace autopas