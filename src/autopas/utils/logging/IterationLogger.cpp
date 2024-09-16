/**
 * @file IterationLogger.cpp
 * @author F. Gratl
 * @date 25/01/2021
 */

#include "IterationLogger.h"

#include "utils/Timer.h"

autopas::IterationLogger::IterationLogger(const std::string &outputSuffix, bool energyMeasurements)
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
  std::string csvHeader =
      "Date,"
      "Iteration,"
      "inTuningPhase,"
      "{},"
      "iteratePairwise[ns],"
      "remainderTraversal[ns],"
      "rebuildNeighborLists[ns],"
      "iteratePairwiseTotal[ns],"
      "tuning[ns]";
  if (energyMeasurements) {
    csvHeader.append(
        ",energyPsys[J],"
        "energyPkg[J],"
        "energyRam[J],"
        "energyGPU[J],"
        "energyCores[J]");
  }
  headerLogger->info(csvHeader, Configuration().getCSVHeader());
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
                                            bool inTuningPhase, long timeTuning,
                                            const IterationMeasurements &measurements) {
#ifdef AUTOPAS_LOG_ITERATIONS
  const auto &[timeIteratePairwise, timeRemainderTraversal, timeRebuild, timeTotal, energyMeasurementsPossible,
               energyPsys, energyPkg, energyRam, energyGPU, energyCores, energyTotal] = measurements;
  if (energyMeasurementsPossible) {
    spdlog::get(_loggerName)
        ->info("{},{},{},{},{},{},{},{},{},{},{},{},{}", iteration, inTuningPhase ? "true" : "false",
               configuration.getCSVLine(), timeIteratePairwise, timeRemainderTraversal, timeRebuild, timeTotal,
               timeTuning, energyPsys, energyPkg, energyRam, energyGPU, energyCores);
  } else {
    spdlog::get(_loggerName)
        ->info("{},{},{},{},{},{},{},{}", iteration, inTuningPhase ? "true" : "false", configuration.getCSVLine(),
               timeIteratePairwise, timeRemainderTraversal, timeRebuild, timeTotal, timeTuning);
  }
#endif
}
