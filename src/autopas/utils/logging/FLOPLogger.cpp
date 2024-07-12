/**
 * @file FLOPLogger.cpp
 * @author F. Gratl & S. Newcome
 * @date 17/04/2024
 */

#include "FLOPLogger.h"

#include <string>

#include "Logger.h"
#include "autopas/utils/Timer.h"

autopas::FLOPLogger::FLOPLogger(const std::string &outputSuffix) : _loggerName("FLOPLogger" + outputSuffix) {
#ifdef AUTOPAS_LOG_FLOPS
  const auto *fillerAfterSuffix = outputSuffix.empty() or outputSuffix.back() == '_' ? "" : "_";
  const auto outputFileName("AutoPas_FLOPCount_" + outputSuffix + fillerAfterSuffix +
                            autopas::utils::Timer::getDateStamp() + ".csv");
  // Start of workaround: Because we want to use an asynchronous logger we can't quickly switch patterns for the header.
  // Create and register a non-asychronous logger to write the header.
  const auto headerLoggerName = _loggerName + "header";
  auto headerLogger = spdlog::basic_logger_mt(headerLoggerName, outputFileName);
  // set the pattern to the message only
  headerLogger->set_pattern("%v");
  // print csv header
  headerLogger->info(
      "Date,"
      "Iteration,"
      "Number FLOPs,"
      "Hit Rate");
  spdlog::drop(headerLoggerName);
  // End of workaround

  // create and register the actual logger
  auto logger = spdlog::basic_logger_mt<spdlog::async_factory>(_loggerName, outputFileName);
  // set pattern to provide date
  logger->set_pattern("%Y-%m-%d %T,%v");
#endif
}

autopas::FLOPLogger::~FLOPLogger() {
#ifdef AUTOPAS_LOG_FLOPS
  spdlog::drop(_loggerName);
#endif
}

void autopas::FLOPLogger::logIteration(size_t iteration, size_t numFLOPs, double hitRate) {
#ifdef AUTOPAS_LOG_FLOPS
  // Convert numeric_limits<size_t>::max and numeric_limits<double>::quiet_NaN to empty strings to represent no implementation
  const auto numFLOPsStr = numFLOPs != std::numeric_limits<size_t>::max() ? std::to_string(numFLOPs) : "Not Implemented";
  const auto hitRateStr = not std::isnan(hitRate) ? std::to_string(hitRate) : "Not Implemented";
  spdlog::get(_loggerName)->info("{},{},{}", iteration, numFLOPsStr, hitRateStr);
#endif
}
