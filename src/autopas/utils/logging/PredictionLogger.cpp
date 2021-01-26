/**
 * @file PredictionLogger.cpp.cc
 * @author F. Gratl
 * @date 25/01/2021
 */

#include "PredictionLogger.h"
autopas::PredictionLogger::PredictionLogger() {
#ifdef AUTOPAS_Log_Predictions
  auto outputFileName("predictions_" + utils::Timer::getDateStamp() + ".csv");
  // create and register the logger
  spdlog::basic_logger_mt<spdlog::async_factory>(loggerName(), outputFileName);
#endif
}
autopas::PredictionLogger::~PredictionLogger() { spdlog::drop(loggerName()); }
void autopas::PredictionLogger::logAllPredictions(
    const std::set<Configuration> &configurations,
    const std::unordered_map<Configuration, size_t, ConfigHash> &configurationPredictions, size_t predictionErrorValue,
    size_t tuningPhaseCounter) {
#ifdef AUTOPAS_Log_Predictions
  spdlog::get(loggerName())->info("Predictions for tuning phase {}", tuningPhaseCounter);
  for (const auto &configuration : configurations) {
    auto prediction = configurationPredictions.at(configuration);
    spdlog::get(loggerName())
        ->info("Prediction for {} : {}", configuration.toString(),
               prediction == predictionErrorValue ? std::to_string(prediction) : "none");
  }
#endif
}
