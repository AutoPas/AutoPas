/**
 * @file IterationLogger.cpp
 * @author F. Gratl
 * @date 25/01/2021
 */

#include "IterationLogger.h"

#include "utils/Timer.h"

autopas::IterationLogger::IterationLogger() {
#ifdef AUTOPAS_Log_Iterations
  auto outputFileName("AutoPas_iterationPerformance_" + utils::Timer::getDateStamp() + ".csv");
  // Start of workaround: Because we want to use an asynchronous logger we can't quickly switch patterns for the header.
  // create and register a non-asychronous logger to write the header
  auto headerLoggerName = loggerName() + "header";
  auto headerLogger = spdlog::basic_logger_mt(headerLoggerName, outputFileName);
  // set the pattern to the message only
  headerLogger->set_pattern("%v");
  // print csv header
  headerLogger->info(
      "Date,Iteration,inTuningPhase,{},iteratePairwise[ns],rebuildNeighborLists[ns],wholeIteration[ns],tuning[ns]",
      Configuration().csvHeader());
  spdlog::drop(headerLoggerName);
  // End of workaround

  // create and register the actual logger
  auto logger = spdlog::basic_logger_mt<spdlog::async_factory>(loggerName(), outputFileName);
  // set pattern to provide date
  logger->set_pattern("%Y-%m-%d %T,%v");
#endif
}

autopas::IterationLogger::~IterationLogger() {
#ifdef AUTOPAS_Log_Iterations
  spdlog::drop(loggerName());
#endif
}

void autopas::IterationLogger::logTimeTuning(long timeTuning) {
#ifdef AUTOPAS_Log_Iterations
  bufferTimeTuning = timeTuning;
#endif
}

void autopas::IterationLogger::logIteration(const autopas::Configuration &configuration, size_t iteration,
                                            bool inTuningPhase, long timeIteratePairwise, long timeRebuildNeighborLists,
                                            long timeWholeIteration) {
#ifdef AUTOPAS_Log_Iterations
  spdlog::get(loggerName())
      ->info("{},{},{},{},{},{},{}", iteration, inTuningPhase ? "true" : "false", configuration.csvLine(),
             timeIteratePairwise, timeRebuildNeighborLists, timeWholeIteration, bufferTimeTuning);

  // reset buffer
  bufferTimeTuning = 0;
#endif
}
