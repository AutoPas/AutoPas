/**
 * @file PredictionLogger.cpp
 * @author F. Gratl
 * @date 25/01/2021
 */

#include "PredictionLogger.h"

#include "autopas/utils/Timer.h"

autopas::PredictionLogger::PredictionLogger(const std::string &outputSuffix)
    : _loggerName("PredictionLogger" + outputSuffix) {
#ifdef AUTOPAS_LOG_PREDICTIONS
  auto outputFileName("AutoPas_predictions_" + outputSuffix + utils::Timer::getDateStamp() + ".csv");
  // Start of workaround: Because we want to use an asynchronous logger we can't quickly switch patterns for the header.
  // create and register a non-asychronous logger to write the header
  auto headerLoggerName = _loggerName + "header";
  auto headerLogger = spdlog::basic_logger_mt(headerLoggerName, outputFileName);
  // set the pattern to the message only
  headerLogger->set_pattern("%v");
  // print csv header
  headerLogger->info("Date,Tuning Phase,{},Prediction", Configuration().getCSVHeader());
  spdlog::drop(headerLoggerName);
  // End of workaround

  // create and register the actual logger
  auto logger = spdlog::basic_logger_mt<spdlog::async_factory>(_loggerName, outputFileName);
  // set pattern to provide date
  logger->set_pattern("%Y-%m-%d %T,%v");
#endif
}

autopas::PredictionLogger::~PredictionLogger() {
#ifdef AUTOPAS_LOG_PREDICTIONS
  spdlog::drop(_loggerName);
#endif
}

void autopas::PredictionLogger::logAllPredictions(
    const std::set<Configuration> &configurations,
    const std::unordered_map<Configuration, size_t, ConfigHash> &configurationPredictions, size_t predictionErrorValue,
    size_t tuningPhaseCounter) {
#ifdef AUTOPAS_LOG_PREDICTIONS
  for (const auto &configuration : configurations) {
    auto prediction = configurationPredictions.at(configuration);
    spdlog::get(_loggerName)
        ->info("{},{},{}", tuningPhaseCounter, configuration.getCSVLine(),
               prediction == predictionErrorValue ? std::to_string(prediction) : "none");
  }
#endif
}
