/**
 * @file Logger.h
 * @author seckler
 * @date 20.04.18
 */

#include <spdlog/sinks/ostream_sink.h>
#include <spdlog/spdlog.h>
#include <iostream>

#pragma once

// this is thread safe
#define AutoPasLogger spdlog::get("AutoPasLog")

namespace autopas {

namespace logger {

static bool created = false;

static void create(std::string& filename) {
  if (created) {
    return;
  }
  spdlog::basic_logger_mt("AutoPasLog", filename);
  created = true;
}

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
