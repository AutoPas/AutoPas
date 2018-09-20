/**
 * @file Logger.h
 * @author seckler
 * @date 20.04.18
 */

#pragma once

#include <spdlog/sinks/ostream_sink.h>
#include <spdlog/spdlog.h>
#include <iostream>

/**
 * this gives you the logger for autopas. call this once the logger has been
 * initialized.
 */
#define AutoPasLogger spdlog::get("AutoPasLog")

/**
 * Returns the filename without full path.
 */
#define __FILENAME__                           \
  ({                                           \
    const char* pStr = strrchr(__FILE__, '/'); \
    pStr ? pStr + 1 : __FILE__;                \
  })

#ifdef AUTOPAS_VERBOSE_LOG
/**
 * Macro for logging providing common meta information.
 * @param lvl Possible levels: trace, debug, info, warn, error, critical.
 * @param fmt Message with formatting tokens
 * @param ... Formatting arguments
 */
#define AutoPasLog(lvl, fmt, ...)                                        \
  {                                                                      \
    size_t textwidth = 26; /* If filenames get cropped increase this! */ \
    std::string s;                                                       \
    s.reserve(textwidth);                                                \
    s.append(__FILENAME__);                                              \
    s.append(":");                                                       \
    s.append(std::to_string(__LINE__));                                  \
    s.resize(textwidth, ' ');                                            \
    AutoPasLogger->lvl("[{}] " fmt, s, ##__VA_ARGS__);                   \
  }
#else
/**
 * Macro for logging providing common meta information without filename.
 * @param lvl Possible levels: trace, debug, info, warn, error, critical.
 * @param fmt Message with formatting tokens
 * @param ... Formatting arguments
 */
#define AutoPasLog(lvl, fmt, ...) \
  { AutoPasLogger->lvl(fmt, ##__VA_ARGS__); }

#endif

namespace autopas {
/**
 * Logger class to provide interface to basic functions of the logger.
 * You can create the spdlog's logger or delete it using the provided functions.
 */
class Logger {
 public:
  /**
   * create a logger writing to the file system
   * @param filename
   */
  static void create(std::string& filename) {
    // drop an already registered Logger if it exists
    if (spdlog::get("AutoPasLog")) spdlog::drop("AutoPasLog");
    spdlog::basic_logger_mt("AutoPasLog", filename);
  }

  /**
   * create a logger with an arbitrary ostream.
   * default is std::cout
   * @param oss
   */
  static void create(std::ostream& oss = std::cout) {
    // drop an already registered Logger if it exists
    if (spdlog::get("AutoPasLog")) spdlog::drop("AutoPasLog");
    auto ostream_sink = std::make_shared<spdlog::sinks::ostream_sink_mt>(oss);
    auto logger = std::make_shared<spdlog::logger>("AutoPasLog", ostream_sink);
    spdlog::register_logger(logger);
  }

  /**
   * removes the logger. This should only be done at teardown of the simulation or
   * for tests.
   * logging after the logger has been removed and no new logger has been defined
   * will lead to undefined behavior!!!
   */
  static void unregister() { spdlog::drop("AutoPasLog"); }

  /**
   * disable the logger
   */
  static void disable() { spdlog::set_level(spdlog::level::off); }
};  // class Logger
}  // namespace autopas
