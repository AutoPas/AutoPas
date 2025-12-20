//
// Created by manish on 12/11/25.
//

#include "GlobalVariableLogger.h"

#include <string>

#include "autopas/utils/Timer.h"
#include "autopas/utils/logging/Logger.h"

GlobalVariableLogger::GlobalVariableLogger(const std::string &outputSuffix)
    : _loggerName("GlobalsLogger" + outputSuffix) {
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
      "Date, "
      "Iteration, "
      "Average Virial, "
      "Virial Sum, "
      "Potential Energy ");
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

void GlobalVariableLogger::logGlobals(const size_t iteration, const double virial, const double potentialEnergy) {
#ifdef MD_FLEXIBLE_CALC_GLOBALS
double avg_virial = virial/iteration;
  spdlog::get(_loggerName)->info("{},{},{},{}", iteration, avg_virial , virial, potentialEnergy);
#endif
}