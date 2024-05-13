/**
 * @file Logger.cpp
 * @author F. Gratl
 * @date 13.05.24
 */
#include "Logger.h"

#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/ostream_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

namespace autopas::Logger {

const std::string &loggerName() {
  // only instantiate the string once
  static const std::string loggerName = "AutoPasLog";
  return loggerName;
}

void create(std::string &filename) {
  // drop an already registered Logger if it exists
  if (spdlog::get(loggerName())) {
    unregister();
  }
  spdlog::basic_logger_mt(loggerName(), filename);
}

void create(std::ostream &oss) {
  // drop an already registered Logger if it exists
  if (spdlog::get(loggerName())) {
    unregister();
  }
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
  spdlog::register_logger(logger);
}

void unregister() { spdlog::drop(loggerName()); }

std::shared_ptr<spdlog::logger> get() { return spdlog::get(loggerName()); }
}  // namespace autopas::Logger