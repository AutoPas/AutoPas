/**
 * @file TuningDataLogger.cpp
 * @author F. Gratl
 * @date 25/01/2021
 */

#include "TuningDataLogger.h"

#include <spdlog/async.h>

#include "utils/Timer.h"

autopas::TuningDataLogger::TuningDataLogger(std::size_t numSamples, const std::string &outputSuffix)
    : _loggerName("TuningDataLogger" + outputSuffix) {
#ifdef AUTOPAS_LOG_TUNINGDATA
  const auto *fillerAfterSuffix = outputSuffix.empty() or outputSuffix.back() == '_' ? "" : "_";
  auto outputFileName("AutoPas_tuningData_" + outputSuffix + fillerAfterSuffix + utils::Timer::getDateStamp() + ".csv");
  // Start of workaround: Because we want to use an asynchronous logger we can't quickly switch patterns for the header.
  // create and register a non-asychronous logger to write the header
  auto headerLoggerName = _loggerName + "header";
  auto headerLogger = spdlog::basic_logger_mt(headerLoggerName, outputFileName);
  // set the pattern to the message only
  headerLogger->set_pattern("%v");
  std::stringstream samplesHeader;
  for (std::size_t i = 0; i < numSamples; ++i) {
    samplesHeader << ",sample" << i;
  }
  // print csv header
  headerLogger->info("Date,Iteration,{}{},Reduced,Smoothed", Configuration().getCSVHeader(), samplesHeader.str());
  spdlog::drop(headerLoggerName);
  // End of workaround

  // create and register the actual logger
  auto logger = spdlog::basic_logger_mt<spdlog::async_factory>(_loggerName, outputFileName);
  // set pattern to provide date
  logger->set_pattern("%Y-%m-%d %T,%v");
#endif
}

autopas::TuningDataLogger::~TuningDataLogger() {
#ifdef AUTOPAS_LOG_TUNINGDATA
  spdlog::drop(_loggerName);
#endif
}

void autopas::TuningDataLogger::logTuningData(const autopas::Configuration &configuration,
                                              const std::vector<long> &samplesRebuildingNeighborLists,
                                              const std::vector<long> &samplesNotRebuildingNeighborLists,
                                              std::size_t iteration, long reducedValue, long smoothedValue) {
#ifdef AUTOPAS_LOG_TUNINGDATA
  spdlog::get(_loggerName)
      ->info("{},{},{},{},{},{}", iteration, configuration.getCSVLine(),
             utils::ArrayUtils::to_string(samplesRebuildingNeighborLists, ",", {"", ""}),
             utils::ArrayUtils::to_string(samplesNotRebuildingNeighborLists, ",", {"", ""}), reducedValue,
             smoothedValue);
#endif
}
