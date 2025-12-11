//
// Created by manish on 12/11/25.
//

#include "GlobalVariableLogger.h"

#include <string>

#include "autopas/utils/logging/Logger.h"
#include "autopas/utils/Timer.h"

GlobalVariableLogger::GlobalVariableLogger(const std::string &outputSuffix) : _loggerName("GlobalsLogger" + outputSuffix) {
#ifdef MD_FLEXIBLE_CALC_GLOBALS
  const auto *fillerAfterSuffix = outputSuffix.empty() or outputSuffix.back() == '_' ? "" : "_";
  const auto outputFileName("MD_Flexible_Globals_" + outputSuffix + fillerAfterSuffix +
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
      "Iteration,"//
      //"Number FLOPs,"
      //"Hit Rate"
      );
  spdlog::drop(headerLoggerName);
  // End of workaround

  // create and register the actual logger
  auto logger = spdlog::basic_logger_mt<spdlog::async_factory>(_loggerName, outputFileName);
  // set pattern to provide date
  logger->set_pattern("%Y-%m-%d %T,%v");
#endif
}

GlobalVariableLogger::~GlobalVariableLogger() {
#ifdef MD_FLEXIBLE_CALC_GLOBALS
  spdlog::drop(_loggerName);
#endif
}

void GlobalVariableLogger::logGlobals(size_t iteration) {
#ifdef MD_FLEXIBLE_CALC_GLOBALS
  spdlog::get(_loggerName)->info("{}", iteration);
#endif
}