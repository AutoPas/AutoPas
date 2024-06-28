/**
 * @file IterationLogger.cpp
 * @author F. Gratl
 * @date 25/01/2021
 */

#include "IterationLogger.h"

#include "utils/Timer.h"

autopas::IterationLogger::IterationLogger(const std::string &outputSuffix)
    : _loggerName("IterationLogger" + outputSuffix) {
#ifdef AUTOPAS_LOG_ITERATIONS
  const auto *fillerAfterSuffix = outputSuffix.empty() or outputSuffix.back() == '_' ? "" : "_";
  const auto outputFileName("AutoPas_iterationPerformance_" + outputSuffix + fillerAfterSuffix +
                            utils::Timer::getDateStamp() + ".csv");
  // Start of workaround: Because we want to use an asynchronous logger we can't quickly switch patterns for the header.
  // Create and register a non-asynchronous logger to write the header.
  const auto headerLoggerName = _loggerName + "header";
  auto headerLogger = spdlog::basic_logger_mt(headerLoggerName, outputFileName);
  // set the pattern to the message only
  headerLogger->set_pattern("%v");
  // print csv header
  headerLogger->info(
      "Date,"
      "Iteration,"
      "Functor,"
      "inTuningPhase,"
      "{},"
      "iterateInteractions[ns],"
      "remainderTraversal[ns],"
      "rebuildNeighborLists[ns],"
      "iterateInteractionsTotal[ns],"
      "tuning[ns],"
      "energyPsys[J],"
      "energyPkg[J],"
      "energyRam[J]",
      Configuration().getCSVHeader());
  spdlog::drop(headerLoggerName);
  // End of workaround

  // create and register the actual logger
  auto logger = spdlog::basic_logger_mt<spdlog::async_factory>(_loggerName, outputFileName);
  // set pattern to provide date
  logger->set_pattern("%Y-%m-%d %T,%v");
#endif
}

autopas::IterationLogger::~IterationLogger() {
#ifdef AUTOPAS_LOG_ITERATIONS
  spdlog::drop(_loggerName);
#endif
}

void autopas::IterationLogger::logIteration(const autopas::Configuration &configuration, size_t iteration,
                                            std::string functorName, bool inTuningPhase, long timeIterateInteractions,
                                            long timeRemainderTraversal, long timeRebuildNeighborLists,
                                            long timeIterateInteractionsTotal, long timeTuning, double energyPsys,
                                            double energyPkg, double energyRam) {
#ifdef AUTOPAS_LOG_ITERATIONS
  spdlog::get(_loggerName)
      ->info("{},{},{},{},{},{},{},{},{},{},{},{}", iteration, functorName, inTuningPhase ? "true" : "false",
             configuration.getCSVLine(), timeIterateInteractions, timeRemainderTraversal, timeRebuildNeighborLists,
             timeIterateInteractionsTotal, timeTuning, energyPsys, energyPkg, energyRam);
#endif
}
