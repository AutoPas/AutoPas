/**
 * @file Logger.h
 * @author seckler
 * @date 20.04.18
 */

#pragma once

#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/ostream_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <iostream>

/**
 * Macro for logging providing common meta information without filename.
 * @param lvl Possible levels: trace, debug, info, warn, error, critical.
 * @param fmt Message with formatting tokens
 * @param ... Formatting arguments
 * @note A ';' is enforced at the end of the macro.
 */
// #define AutoPasLog(lvl, fmt, ...) SPDLOG_LOGGER_##lvl(spdlog::get("AutoPasLog"), fmt, ##__VA_ARGS__)

#define AutoPasLog(lvl, fmt, ...) \
do { \
auto _logger_ptr = spdlog::get("AutoPasLog"); \
if (!_logger_ptr) { \
fprintf(stderr, "CRASH PREVENTED: Null logger in %s:%d\n", __FILE__, __LINE__); \
} else { \
SPDLOG_LOGGER_##lvl(_logger_ptr, fmt, ##__VA_ARGS__); \
} \
} while (0)

namespace autopas {
/**
 * Logger class to provide interface to basic functions of the logger.
 * You can create the spdlog's logger or delete it using the provided functions.
 */
class Logger {
  static inline auto loggerName() { return "AutoPasLog"; };

 public:
  /**
   * Constructor for a Logger that logs to an output stream. Per default to std::cout.
   * @param logOutputStream
   */
  explicit Logger(std::ostream &logOutputStream = std::cout) {
    // Only first AutoPas instance creates the logger
    if (_instanceCount.fetch_add(1) == 0) {
      create(logOutputStream);
    }
  }

  /**
   * Constructor for a Logger that logs to a file.
   * @param filename
   */
  explicit Logger(const std::string &filename) {
    // Only first AutoPas instance creates the logger
    if (_instanceCount.fetch_add(1) == 0) {
      create(filename);
    }
  }

  /**
   * Destructor that shuts down the logger if this was the last instance.
   */
  ~Logger() {
    // Only the last AutoPas instance shuts down the logger
    if (_instanceCount.fetch_sub(1) == 1) {
      unregister();
    }
  }

  /**
   * Typalias for log levels.
   */
  using LogLevel = spdlog::level::level_enum;

  /**
   * Get a pointer to the actual logger object.
   * @return Pointer to logger.
   */
  static auto get() { return spdlog::get(loggerName()); }

 private:
  /**
   * Atomic instance counter to avoid shutting down the AutoPas logger when another instance is still running.
   */
  static inline std::atomic<size_t> _instanceCount{0};

  /**
   * Create a logger with an arbitrary ostream.
   * Default is std::cout
   * @param oss
   */
  static void create(std::ostream &oss = std::cout) {
    // drops an already registered Logger if it exists
    unregister();

    std::shared_ptr<spdlog::sinks::sink> ostreamSink;
#ifdef AUTOPAS_COLORED_CONSOLE_LOGGING
    if (oss.rdbuf() == std::cout.rdbuf()) {
      ostreamSink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    } else if (oss.rdbuf() == std::cerr.rdbuf()) {
      ostreamSink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
    } else {  // no color for streams other than cout / cerr
      ostreamSink = std::make_shared<spdlog::sinks::ostream_sink_mt>(oss);
    }
#else
    ostreamSink = std::make_shared<spdlog::sinks::ostream_sink_mt>(oss);
#endif
    auto logger = std::make_shared<spdlog::logger>(loggerName(), ostreamSink);
    // The logger is normally only flushed on successful program termination.
    // This line ensures flushing when log messages of level warning or more severe are created.
    logger->flush_on(spdlog::level::warn);
    spdlog::register_logger(logger);
  }

  /**
   * Static: Create a file logger
   */
  static void create(const std::string &filename) {
    unregister();
    // basic_logger_mt creates a simple file sink
    auto logger = spdlog::basic_logger_mt(loggerName(), filename);
    logger->flush_on(spdlog::level::warn);
  }

  /**
   * Removes the logger. This should only be done at teardown of the simulation or
   * for tests.
   * Logging after the logger has been removed and no new logger has been defined
   * will lead to undefined behavior!
   */
  static void unregister() { spdlog::drop(loggerName()); }
};  // class Logger
}  // namespace autopas
