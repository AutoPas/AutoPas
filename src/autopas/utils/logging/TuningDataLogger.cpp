/**
 * @file TuningDataLogger.cpp
 * @author F. Gratl
 * @date 25/01/2021
 */

#include "TuningDataLogger.h"

#include "utils/Timer.h"

autopas::TuningDataLogger::TuningDataLogger(size_t numSamples, const std::string &outputSuffix) {
#ifdef AUTOPAS_LOG_TUNINGDATA
  auto outputFileName("AutoPas_tuningData_" + outputSuffix + utils::Timer::getDateStamp() + ".csv");
  // Start of workaround: Because we want to use an asynchronous logger we can't quickly switch patterns for the header.
  // create and register a non-asychronous logger to write the header
  auto headerLoggerName = loggerName() + "header";
  auto headerLogger = spdlog::basic_logger_mt(headerLoggerName, outputFileName);
  // set the pattern to the message only
  headerLogger->set_pattern("%v");
  std::stringstream samplesHeader;
  for (int i = 0; i < numSamples; ++i) {
    samplesHeader << ",sample" << i;
  }
  // print csv header
  headerLogger->info("Date,Iteration,{}{},Reduced,Smoothed", Configuration().getCSVHeader(), samplesHeader.str());
  spdlog::drop(headerLoggerName);
  // End of workaround

  // create and register the actual logger
  auto logger = spdlog::basic_logger_mt<spdlog::async_factory>(loggerName(), outputFileName);
  // set pattern to provide date
  logger->set_pattern("%Y-%m-%d %T,%v");
#endif
}

autopas::TuningDataLogger::~TuningDataLogger() {
#ifdef AUTOPAS_LOG_TUNINGDATA
  spdlog::drop(loggerName());
#endif
}

void autopas::TuningDataLogger::logTuningData(const autopas::Configuration &configuration,
                                              const std::vector<size_t> &samples, size_t iteration, size_t reducedValue,
                                              size_t smoothedVale) {
#ifdef AUTOPAS_LOG_TUNINGDATA
  spdlog::get(loggerName())
      ->info("{},{},{},{},{}", iteration, configuration.getCSVLine(),
             utils::ArrayUtils::to_string(samples, ",", {"", ""}), reducedValue, smoothedVale);
#endif
}
