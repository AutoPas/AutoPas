/**
 * @file FLOPLogger.h
 * @author F. Gratl & S. Newcome
 * @date 17/04/2024
 */

#pragma once

#include <spdlog/async.h>
#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/ostream_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace autopas {

/**
 * Helper to log FLOP count and HitRate for AutoPas::iteratePairwise() calls with the functors in the molecular dynamics
 * library to a csv file for easier analysis.
 *
 * It uses an asynchronous spd logger to write a csv file named "AutoPas_FLOPCount_<dateStamp>.csv".
 *
 * By default logging the data is disabled. It can be enabled by setting the cmake variable MD_FLEXIBLE_LOG_FLOPS
 * to ON.
 */
class FLOPLogger {
 public:
  /**
   * Constructor initializes the logger and sets the output file name.
   * @param outputSuffix Suffix for all output files produced by this class.
   */
  explicit FLOPLogger(const std::string &outputSuffix = "");

  /**
   * Destructor drops the logger from the spd registry.
   */
  ~FLOPLogger();

  /**
   * Log the given arguments and the internal buffer to the csv file.
   * @param iteration
   * @param numFLOPs number of FLOPs
   * @param hitRate percentage of distance calculations that result in force contributions.
   */
  void logIteration(size_t iteration, size_t numFLOPs, double hitRate);

 private:
  std::string _loggerName;
};

}  // namespace autopas
