/**
 * @file IterationLogger.cpp
 * @author F. Gratl
 * @date 25/01/2021
 */

#include "IterationLogger.h"

#include "utils/Timer.h"

autopas::IterationLogger::IterationLogger(const std::string &outputSuffix, bool energyMeasurements,
                                          bool usingFastParticlesBuffer)
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
      "Functor,"
      "inTuningPhase,"
      "{},"
      "computeInteractions[ns],"
      "remainderTraversal[ns],"
      "rebuildNeighborLists[ns],"
      "computeInteractionsTotal[ns],"
      "tuning[ns]";
  if (usingFastParticlesBuffer) {
    csvHeader.append(
        ",numParticlesBuffer,"
        "numParticlesOwned,"
        "numParticlesHalo,"
        "numParticlesFast,"
        "estimateNumParticlesBuffer,"
        "estimateRemainderTraversal[ns],"
        "estimateRebuildNeighborLists[ns],"
        "doDynamicRebuild");
  }
  if (energyMeasurements) {
    csvHeader.append(
        ",energyWatts[W],"
        "energyJoules[J],"
        "energyDeltaT[s]");
  }
  headerLogger->info(fmt::runtime(csvHeader), Configuration().getCSVHeader());
  spdlog::drop(headerLoggerName);
  // End of workaround

  // create and register the actual logger
  auto logger = spdlog::basic_logger_mt<spdlog::async_factory>(_loggerName, outputFileName);
  // set pattern to provide date
  logger->set_pattern("%Y-%m-%d %T,%v");
#endif
}

void autopas::IterationLogger::logIteration(const autopas::Configuration &configuration, size_t iteration,
                                            const std::string &functorName, bool inTuningPhase, long timeTuning,
                                            const IterationMeasurements &measurements,
                                            bool usingFastParticlesBuffer) const {
#ifdef AUTOPAS_LOG_ITERATIONS
  const auto &[timeIteratePairwise, timeRemainderTraversal, timeRebuild, timeTotal, energyMeasurementsPossible,
               energyWatts, energyJoules, energyDeltaT, energyTotalRebuild, energyTotalNonRebuild, energyTotal,
               numParticlesBuffer, numParticlesOwned, numParticlesHalo, numParticlesFast, estimateNumParticlesBuffer,
               estimateRemainderTraversalTime, estimateRebuildNeighborTime, doDynamicRebuild] = measurements;
  fmt::memory_buffer logBuffer;

  fmt::format_to(std::back_inserter(logBuffer), "{},{},{},{},{},{},{},{},{}", iteration, functorName,
                 inTuningPhase ? "true" : "false", configuration.getCSVLine(), timeIteratePairwise,
                 timeRemainderTraversal, timeRebuild, timeTotal, timeTuning);

  if (usingFastParticlesBuffer) {
    fmt::format_to(std::back_inserter(logBuffer), ",{},{},{},{},{},{},{},{}", numParticlesBuffer, numParticlesOwned,
                   numParticlesHalo, numParticlesFast, estimateNumParticlesBuffer, estimateRemainderTraversalTime,
                   estimateRebuildNeighborTime, doDynamicRebuild ? "true" : "false");
  }

  if (energyMeasurementsPossible) {
    fmt::format_to(std::back_inserter(logBuffer), ",{},{},{}", energyWatts, energyJoules, energyDeltaT);
  }

  spdlog::get(_loggerName)->info("{}", fmt::to_string(logBuffer));
#endif
}

autopas::IterationLogger::~IterationLogger() {
#ifdef AUTOPAS_LOG_ITERATIONS
  spdlog::drop(_loggerName);
#endif
}