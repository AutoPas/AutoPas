/**
 * @file TuningResultLogger.cpp
 * @author F. Gratl
 * @date 25/01/2021
 */

#include "TuningResultLogger.h"

#include "utils/Timer.h"

autopas::TuningResultLogger::TuningResultLogger() {
#ifdef AUTOPAS_Log_TuningResults
  auto outputFileName("AutoPas_tuningResults_" + utils::Timer::getDateStamp() + ".csv");
  // Start of workaround: Because we want to use an asynchronous logger we can't quickly switch patterns for the header.
  // create and register a non-asychronous logger to write the header
  auto headerLoggerName = loggerName() + "header";
  auto headerLogger = spdlog::basic_logger_mt(headerLoggerName, outputFileName);
  // set the pattern to the message only
  headerLogger->set_pattern("%v");
  // print csv header
  headerLogger->info("Date,Iteration,{},tuning[ns]", Configuration().csvHeader());
  spdlog::drop(headerLoggerName);
  // End of workaround

  // create and register the actual logger
  auto logger = spdlog::basic_logger_mt<spdlog::async_factory>(loggerName(), outputFileName);
  // set pattern to provide date
  logger->set_pattern("%Y-%m-%d %T,%v");
#endif
}

autopas::TuningResultLogger::~TuningResultLogger() {
#ifdef AUTOPAS_Log_TuningResults
  spdlog::drop(loggerName());
#endif
}

void autopas::TuningResultLogger::logTuningResult(const autopas::Configuration &configuration, size_t iteration,
                                                  long timeTuning) {
#ifdef AUTOPAS_Log_TuningResults
  spdlog::get(loggerName())->info("{},{},{}", iteration, configuration.csvLine(), timeTuning);
#endif
}
