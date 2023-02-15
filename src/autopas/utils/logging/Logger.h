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

#ifdef AUTOPAS_VERBOSE_LOG
/**
 * Helper Macro to get only the basename
 */
#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

/**
 * Macro for logging providing common meta information.
 * @param lvl Possible levels: trace, debug, info, warn, error, critical.
 * @param fmt Message with formatting tokens
 * @param ... Formatting arguments
 */
#define AutoPasLog(lvl, fmt, ...)                                        \
  {                                                                      \
    size_t textwidth = 35; /* If filenames get cropped increase this! */ \
    std::string s;                                                       \
    s.reserve(textwidth);                                                \
    s.append(__FILENAME__);                                              \
    s.append(":");                                                       \
    s.append(std::to_string(__LINE__));                                  \
    s.resize(textwidth, ' ');                                            \
    spdlog::get("AutoPasLog")->lvl("[{}] " fmt, s, ##__VA_ARGS__);       \
  }
#else
/**
 * Macro for logging providing common meta information without filename.
 * @param lvl Possible levels: trace, debug, info, warn, error, critical.
 * @param fmt Message with formatting tokens
 * @param ... Formatting arguments
 * @note A ';' is enforced at the end of the macro.
 */
#define AutoPasLog(lvl, fmt, ...) SPDLOG_LOGGER_##lvl(spdlog::get("AutoPasLog"), fmt, ##__VA_ARGS__)
#endif

namespace autopas {
/**
 * Logger class to provide interface to basic functions of the logger.
 * You can create the spdlog's logger or delete it using the provided functions.
 */
class Logger {
 private:
  static inline auto loggerName() { return "AutoPasLog"; };

 public:
  /**
   * Typalias for log levels.
   */
  using LogLevel = spdlog::level::level_enum;
  /**
   * create a logger writing to the file system
   * @param filename
   */
  static void create(std::string &filename) {
    // drop an already registered Logger if it exists
    if (spdlog::get(loggerName())) spdlog::drop(loggerName());
    spdlog::basic_logger_mt(loggerName(), filename);
  }

  /**
   * create a logger with an arbitrary ostream.
   * default is std::cout
   * @param oss
   */
  static void create(std::ostream &oss = std::cout) {
    // drop an already registered Logger if it exists
    if (spdlog::get(loggerName())) spdlog::drop(loggerName());
    std::shared_ptr<spdlog::sinks::sink> ostream_sink;
#ifdef AUTOPAS_COLORED_CONSOLE_LOGGING
    if (oss.rdbuf() == std::cout.rdbuf()) {
      ostream_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    } else if (oss.rdbuf() == std::cerr.rdbuf()) {
      ostream_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
    } else {  // no color for streams other than cout / cerr
      ostream_sink = std::make_shared<spdlog::sinks::ostream_sink_mt>(oss);
    }
#else
    ostream_sink = std::make_shared<spdlog::sinks::ostream_sink_mt>(oss);
#endif
    auto logger = std::make_shared<spdlog::logger>(loggerName(), ostream_sink);
    spdlog::register_logger(logger);
  }

  /**
   * removes the logger. This should only be done at teardown of the simulation or
   * for tests.
   * logging after the logger has been removed and no new logger has been defined
   * will lead to undefined behavior!!!
   */
  static void unregister() { spdlog::drop(loggerName()); }

  /**
   * Get a pointer to the actual logger object.
   * @return Pointer to logger.
   */
  static auto get() { return spdlog::get(loggerName()); }
};  // class Logger
}  // namespace autopas
