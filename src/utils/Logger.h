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
 * this gives you the logger for autopas. call this once the logger has been initialized.
 */
#define AutoPasLogger spdlog::get("AutoPasLog")

namespace autopas {

namespace logger {

static bool created = false;

/**
 * create a logger writing to the file system
 * @param filename
 */
static void create(std::string& filename) {
  if (created) {
    return;
  }
  spdlog::basic_logger_mt("AutoPasLog", filename);
  created = true;
}

/**
 * create a logger with an arbitrary ostream.
 * default is std::cout
 * @param oss
 */
static void create(std::ostream& oss = std::cout) {
  if (created) {
    return;
  }

  auto ostream_sink = std::make_shared<spdlog::sinks::ostream_sink_mt>(oss);
  auto logger = std::make_shared<spdlog::logger>("AutoPasLog", ostream_sink);
  spdlog::register_logger(logger);
  created = true;
}

/**
 * disable the logger
 */
static void disable() { spdlog::set_level(spdlog::level::off); }
}  // namespace logger
}  // namespace autopas
