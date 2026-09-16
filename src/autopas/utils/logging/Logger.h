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
 * @name Per-level building blocks of AutoPasLog
 * These mirror spdlog's own SPDLOG_LOGGER_<LEVEL> wrappers, but are gated on AUTOPAS_ACTIVE_LEVEL. Levels
 * below the configured one expand to (void)0, so they cost nothing at run time.
 * The purpose of this is to allow for a consumer to keep AutoPas's log level above that of their own spdlog, so they
 * get e.g. their own trace/debug logs without getting ours.
 * @{
 */
#if AUTOPAS_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
#define AutoPasLog_TRACE(...) SPDLOG_LOGGER_CALL(autopas::Logger::get(), spdlog::level::trace, __VA_ARGS__)
#else
#define AutoPasLog_TRACE(...) (void)0
#endif

#if AUTOPAS_ACTIVE_LEVEL <= SPDLOG_LEVEL_DEBUG
#define AutoPasLog_DEBUG(...) SPDLOG_LOGGER_CALL(autopas::Logger::get(), spdlog::level::debug, __VA_ARGS__)
#else
#define AutoPasLog_DEBUG(...) (void)0
#endif

#if AUTOPAS_ACTIVE_LEVEL <= SPDLOG_LEVEL_INFO
#define AutoPasLog_INFO(...) SPDLOG_LOGGER_CALL(autopas::Logger::get(), spdlog::level::info, __VA_ARGS__)
#else
#define AutoPasLog_INFO(...) (void)0
#endif

#if AUTOPAS_ACTIVE_LEVEL <= SPDLOG_LEVEL_WARN
#define AutoPasLog_WARN(...) SPDLOG_LOGGER_CALL(autopas::Logger::get(), spdlog::level::warn, __VA_ARGS__)
#else
#define AutoPasLog_WARN(...) (void)0
#endif

#if AUTOPAS_ACTIVE_LEVEL <= SPDLOG_LEVEL_ERROR
#define AutoPasLog_ERROR(...) SPDLOG_LOGGER_CALL(autopas::Logger::get(), spdlog::level::err, __VA_ARGS__)
#else
#define AutoPasLog_ERROR(...) (void)0
#endif

#if AUTOPAS_ACTIVE_LEVEL <= SPDLOG_LEVEL_CRITICAL
#define AutoPasLog_CRITICAL(...) \
  SPDLOG_LOGGER_CALL(autopas::Logger::get(), spdlog::level::critical, __VA_ARGS__)
#else
#define AutoPasLog_CRITICAL(...) (void)0
#endif
/** @} */

/**
 * Macro for logging providing common meta information without filename.
 * @param lvl Possible levels: TRACE, DEBUG, INFO, WARN, ERROR, CRITICAL.
 * @param fmt Message with formatting tokens
 * @param ... Formatting arguments
 * @note A ';' is enforced at the end of the macro.
 */
#define AutoPasLog(lvl, fmt, ...) AutoPasLog_##lvl(fmt, ##__VA_ARGS__)

namespace autopas {
/**
 * Logger class to provide interface to basic functions of the logger.
 * It manages the lifecycle safely using lazy initialization.
 */
class Logger {
  static inline const std::string loggerName = "AutoPasLog";

 public:
  /**
   * Typealias for log levels.
   */
  using LogLevel = spdlog::level::level_enum;

  /**
   * Get the pointer to the spdlog logger.
   * If the logger does not exist yet (e.g., during static initialization),
   * this function will automatically create a default one (std::cout).
   * This guarantees that Logger::get() never returns nullptr.
   * @return Shared pointer to the logger.
   */
  static std::shared_ptr<spdlog::logger> get() {
    auto logger = spdlog::get(loggerName);
    if (logger) {
      return logger;
    }

    // Logger not found? Create a default one safely.
    std::lock_guard<std::mutex> lock(_registrationMutex);

    // Double-check after locking in case another thread created it just now
    logger = spdlog::get(loggerName);
    if (not logger) {
      logger = create_internal(std::cout);  // Default to cout
    }
    return logger;
  }

  /**
   * Explicitly initialize/reset the logger to write to an output stream.
   * @param logOutputStream Stream to write to (default std::cout)
   */
  static void create(std::ostream &logOutputStream = std::cout) {
    std::lock_guard<std::mutex> lock(_registrationMutex);
    unregister_internal();  // Drop old one first
    create_internal(logOutputStream);
  }

  /**
   * Explicitly initialize/reset the logger to write to a file.
   * @param filename Path to the log file.
   */
  static void create(const std::string &filename) {
    std::lock_guard<std::mutex> lock(_registrationMutex);
    unregister_internal();  // Drop old one first

    auto logger = spdlog::basic_logger_mt(loggerName, filename);
    logger->flush_on(spdlog::level::warn);
    // spdlog automatically registers the logger upon creation via factory methods
  }

  /**
   * Removes the logger from the registry.
   * Usually not needed unless you want to silence leak detectors at the very end of main().
   */
  static void unregister() {
    std::lock_guard<std::mutex> lock(_registrationMutex);
    unregister_internal();
  }

 private:
  /**
   * Mutex to protect creation and destruction of the global logger.
   */
  static inline std::mutex _registrationMutex;

  /**
   * Internal helper to create the logger (must be called with lock held).
   */
  static std::shared_ptr<spdlog::logger> create_internal(std::ostream &oss) {
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
    auto logger = std::make_shared<spdlog::logger>(loggerName, ostreamSink);
    // The logger is normally only flushed on successful program termination.
    // This line ensures flushing when log messages of level warning or more severe are created.
    logger->flush_on(spdlog::level::warn);
    spdlog::register_logger(logger);
    return logger;
  }

  /**
   * Internal helper to drop the logger (must be called with lock held).
   */
  static void unregister_internal() { spdlog::drop(loggerName); }
};  // class Logger
}  // namespace autopas
