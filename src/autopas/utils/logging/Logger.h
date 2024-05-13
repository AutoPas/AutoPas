/**
 * @file Logger.h
 * @author seckler
 * @date 20.04.18
 */

#pragma once

#include <spdlog/spdlog.h>

#include <iostream>

/**
 * Macro for logging providing common meta information without filename.
 * @param lvl Possible levels: trace, debug, info, warn, error, critical.
 * @param fmt Message with formatting tokens
 * @param ... Formatting arguments
 * @note A ';' is enforced at the end of the macro.
 */
#define AutoPasLog(lvl, fmt, ...) SPDLOG_LOGGER_##lvl(spdlog::get("AutoPasLog"), fmt, ##__VA_ARGS__)

namespace autopas {
/**
 * Logger class to provide interface to basic functions of the logger.
 * You can create the spdlog's logger or delete it using the provided functions.
 */
namespace Logger {

/**
 * Returns the static name of the logger.
 * @return
 */
const std::string &loggerName();

/**
 * Type alias for log levels.
 */
using LogLevel = spdlog::level::level_enum;

/**
 * Create a logger writing to the file system.
 * @param filename
 */
void create(std::string &filename);

/**
 * Create a logger with an arbitrary ostream.
 * The default is std::cout.
 * @param oss
 */
void create(std::ostream &oss = std::cout);

/**
 * Removes the logger. This should only be done at teardown of the simulation or for tests.
 * Logging after the logger has been removed and no new logger has been defined will lead to undefined behavior!
 */
void unregister();

/**
 * Get a pointer to the actual logger object.
 * @return Pointer to logger.
 */
std::shared_ptr<spdlog::logger> get();
}  // namespace Logger
}  // namespace autopas
