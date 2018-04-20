/**
 * @file Logger.h
 * @author seckler
 * @date 20.04.18
 */

#include <spdlog/sinks/ostream_sink.h>
#include <spdlog/spdlog.h>
#include <iostream>

#pragma once

/**
 * this gives you the logger for autopas. call this once the logger has been
 * initialized.
 */
#define AutoPasLogger spdlog::get("AutoPasLog")

namespace autopas {

namespace logger {

/**
 * create a logger writing to the file system
 * @param filename
 */
static void create(std::string& filename) {
  spdlog::basic_logger_mt("AutoPasLog", filename);
}

/**
 * create a logger with an arbitrary ostream.
 * default is std::cout
 * @param oss
 */
static void create(std::ostream& oss = std::cout) {
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
}  // namespace logger
}  // namespace autopas
