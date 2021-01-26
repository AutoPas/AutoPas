/**
 * @file PredictionLogger.cpp
 * @author F. Gratl
 * @date 25/01/2021
 */

#include "PredictionLogger.h"

#include "autopas/utils/Timer.h"

autopas::PredictionLogger::PredictionLogger() {
#ifdef AUTOPAS_Log_Predictions
  auto outputFileName("AutoPas_predictions_" + utils::Timer::getDateStamp() + ".csv");
  // create and register the logger
  auto logger = spdlog::basic_logger_mt<spdlog::async_factory>(loggerName(), outputFileName);
  // set the pattern to the message only
  logger->set_pattern("%v");
  logger->info("Tuning phase,{},Prediction", Configuration().csvHeader());
#endif
}

autopas::PredictionLogger::~PredictionLogger() {
#ifdef AUTOPAS_Log_Predictions
  spdlog::drop(loggerName());
#endif
}

void autopas::PredictionLogger::logAllPredictions(
    const std::set<Configuration> &configurations,
    const std::unordered_map<Configuration, size_t, ConfigHash> &configurationPredictions, size_t predictionErrorValue,
    size_t tuningPhaseCounter) {
#ifdef AUTOPAS_Log_Predictions
  for (const auto &configuration : configurations) {
    auto prediction = configurationPredictions.at(configuration);
    spdlog::get(loggerName())
        ->info("{},{},{}", tuningPhaseCounter, configuration.csvLine(),
               prediction == predictionErrorValue ? std::to_string(prediction) : "none");
  }
#endif
}
