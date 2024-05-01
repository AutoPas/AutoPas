/**
 * @file LiveInfoLogger.cpp
 * @author Manuel Lerchner
 * @date 29/04/2024
 */

#include "LiveInfoLogger.h"

#include "utils/Timer.h"

autopas::LiveInfoLogger::LiveInfoLogger(const std::string &outputSuffix)
    : _loggerName("LiveInfoLogger" + outputSuffix) {
#ifdef AUTOPAS_LOG_LIVEINFO
  const auto *fillerAfterSuffix = outputSuffix.empty() or outputSuffix.back() == '_' ? "" : "_";
  _outputFileName =
      "AutoPas_liveInfoLogger_" + outputSuffix + fillerAfterSuffix + utils::Timer::getDateStamp() + ".csv";

  // create and register the actual logger
  auto logger = spdlog::basic_logger_mt<spdlog::async_factory>(_loggerName, _outputFileName);
  // set pattern to provide date
  logger->set_pattern("%Y-%m-%d %T,%v");
  // since this logger only writes rarely flush instantly in order to not lose any information if autopas is killed
  logger->flush_on(spdlog::level::trace);
#endif
}

autopas::LiveInfoLogger::~LiveInfoLogger() {
#ifdef AUTOPAS_LOG_LIVEINFO
  spdlog::drop(_loggerName);
#endif
}

void autopas::LiveInfoLogger::logLiveInfo(const LiveInfo &liveInfo, size_t iteration) {
#ifdef AUTOPAS_LOG_LIVEINFO
  const auto [csvHeader, csvLine] = liveInfo.getCSVLine();

  if (not headerWritten) {
    // Start of workaround: Because we want to use an asynchronous logger we can't quickly switch patterns for the
    // header. Create and register a non-asychronous logger to write the header.
    const auto headerLoggerName = _loggerName + "header";
    auto headerLogger = spdlog::basic_logger_mt(headerLoggerName, _outputFileName);
    // set the pattern to the message only
    headerLogger->set_pattern("%v");
    // print csv header
    headerLogger->info(
        "Date,"
        "Iteration,"
        "{}",
        csvHeader);
    spdlog::drop(headerLoggerName);
    // End of workaround
    headerWritten = true;
  }

  spdlog::get(_loggerName)->info("{},{}", iteration, csvLine);
#endif
}
