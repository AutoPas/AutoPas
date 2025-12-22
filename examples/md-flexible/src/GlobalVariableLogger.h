/**
 * @file GlobalVariableLogger.h
 * @author Manish Mishra
 * @date 12.11.2025
 */

#pragma once

#include <spdlog/async.h>
#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/ostream_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

/**
 * Helper to global variables calculated by the Functor, eg. Virial, Potential Energy, etc.
 *
 * It uses an asynchronous spd logger to write a csv file named "MD_FLEXIBLE_GLOBAL_<dateStamp>.csv".
 *
 * By default logging the global data is disabled. It can be enabled by setting the cmake variable
 * MD_FLEXIBLE_CALC_GLOBALS to ON.
 *
 * When enabled and used with a functor where calculating globals is not implemented (in which case the functor will
 * return the default nonsensical values), "Not Implemented" is outputted instead.
 *
 */
class GlobalVariableLogger {
 public:
  /**
   * Constructor initializes the logger and sets the output file name.
   * @param outputSuffix Suffix for all output files produced by this class.
   */
  explicit GlobalVariableLogger(const std::string &outputSuffix = "");

  /**
   * Destructor drops the logger from the spd registry.
   */
  ~GlobalVariableLogger();

  /**
   * Log the given arguments and the internal buffer to the csv file. If a value is not valid, it is interpreted that
   * the functor has not implemented the relevant function.
   *
   * @param iteration
   * @param virial This includes the virial sum contribution from both 2B and 3B interactions at the current timestep.
   * @param potentialEnergy This includes the potential energy contribution from both 2B and 3B interactions at the
   * current timestep.
   */
  void logGlobals(const size_t iteration, const double virial, const double potentialEnergy);

 private:
  std::string _loggerName;
};
