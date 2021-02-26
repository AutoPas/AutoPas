/**
 * @file TuningResultLogger.cpp
 * @author F. Gratl
 * @date 25/01/2021
 */

#include "TuningResultLogger.h"

#include "utils/Timer.h"

autopas::TuningResultLogger::TuningResultLogger(const std::string &outputSuffix)
    : _loggerName("TuningResultLogger" + outputSuffix) {
#ifdef AUTOPAS_LOG_TUNINGRESULTS
  auto outputFileName("AutoPas_tuningResults_" + outputSuffix + utils::Timer::getDateStamp() + ".csv");
  // create and register the actual logger
  auto logger = spdlog::basic_logger_mt(_loggerName, outputFileName);
  // since this logger only writes rarely flush instantly in order to not lose any information if autopas is killed
  logger->flush_on(spdlog::level::trace);
  // set the pattern to the message only
  logger->set_pattern("%v");
  // print csv header
  logger->info("Date,Iteration,{},tuning[ns]", Configuration().getCSVHeader());
  // set pattern to provide date
  logger->set_pattern("%Y-%m-%d %T,%v");
#endif
}

autopas::TuningResultLogger::~TuningResultLogger() {
#ifdef AUTOPAS_LOG_TUNINGRESULTS
  spdlog::drop(_loggerName);
#endif
}

void autopas::TuningResultLogger::logTuningResult(const autopas::Configuration &configuration, size_t iteration,
                                                  long timeTuning) {
#ifdef AUTOPAS_LOG_TUNINGRESULTS
  spdlog::get(_loggerName)->info("{},{},{}", iteration, configuration.getCSVLine(), timeTuning);
#endif
}
