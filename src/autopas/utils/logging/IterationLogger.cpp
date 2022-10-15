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
  auto outputFileName("AutoPas_iterationPerformance_" + outputSuffix + utils::Timer::getDateStamp() + ".csv");
  // Start of workaround: Because we want to use an asynchronous logger we can't quickly switch patterns for the header.
  // Create and register a non-asychronous logger to write the header.
  auto headerLoggerName = _loggerName + "header";
  auto headerLogger = spdlog::basic_logger_mt(headerLoggerName, outputFileName);
  // set the pattern to the message only
  headerLogger->set_pattern("%v");
  // print csv header
  headerLogger->info(
      "Date,Iteration,inTuningPhase,{},iteratePairwise[ns],rebuildNeighborLists[ns],wholeIteration[ns],tuning[ns],energyPsys,energyPkg,energyRam",
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

void autopas::IterationLogger::logTimeTuning(long timeTuning) {
#ifdef AUTOPAS_LOG_ITERATIONS
  _bufferTimeTuning = timeTuning;
#endif
}

void autopas::IterationLogger::logIteration(const autopas::Configuration &configuration, size_t iteration,
                                            bool inTuningPhase, long timeIteratePairwise, long timeRebuildNeighborLists,
                                            long timeWholeIteration, double energyPsys, double energyPkg,
                                            double energyRam) {
#ifdef AUTOPAS_LOG_ITERATIONS
  spdlog::get(_loggerName)
      ->info("{},{},{},{},{},{},{}", iteration, inTuningPhase ? "true" : "false", configuration.getCSVLine(),
             timeIteratePairwise, timeRebuildNeighborLists, timeWholeIteration, _bufferTimeTuning, energyPsys,
             energyPkg, energyRam);

  // reset buffer
  _bufferTimeTuning = 0;
#endif
}
