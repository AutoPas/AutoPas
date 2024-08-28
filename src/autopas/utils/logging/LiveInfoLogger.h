/**
 * @file LiveInfoLogger.h
 * @author Manuel Lerchner
 * @date 29/04/2024
 */

#pragma once

#include <spdlog/async.h>

#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/tuningStrategy/LiveInfo.h"

namespace autopas {
/**
 * Helper to log the collected LiveInfo data during tuning to a csv file for easier analysis.
 *
 * It uses an asynchronous spd logger to write a csv file named "AutoPas_liveInfoLogger_<dateStamp>.csv".
 *
 * By default logging the data is disabled. It can be enabled by setting the cmake variable AUTOPAS_LOG_LIVEINFO
 * to ON.
 */
class LiveInfoLogger {
 public:
  /**
   * Constructor initializes the logger and sets the output file name.
   * @param outputSuffix Suffix for all output files produced by this class.
   */
  explicit LiveInfoLogger(const std::string &outputSuffix = "");

  /**
   * Destructor drops the logger from the spd registry.
   */
  ~LiveInfoLogger();

  /**
   * Log the given arguments and the internal buffer to the csv file.
   * @param liveInfo The liveInfo Struct to log.
   * @param iteration The iteration number.
   */
  void logLiveInfo(const LiveInfo &liveInfo, size_t iteration);

 private:
  std::string _loggerName;
  std::string _outputFileName;
#ifdef AUTOPAS_LOG_LIVEINFO
  bool headerWritten = false;
#endif
};
}  // namespace autopas
