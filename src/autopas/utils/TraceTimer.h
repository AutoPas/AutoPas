/**
 * @file TraceTimer.h
 * @date 12.02.2026
 * @author muehlhaeusser
 */

#pragma once

#include <spdlog/spdlog.h>

#include "autopas/utils/Timer.h"

namespace autopas::utils {

/**
 * A wrapper around autopas::utils::Timer that only compiles implementation logic
 * if the SPDLOG_ACTIVE_LEVEL is set to TRACE.
 * For a higher log level, functions are empty and will be removed by the compiler.
 */
class TraceTimer {
 public:
  /**
   * @copydoc autopas::utils::Timer::start()
   */
  void start() {
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    _timer.start();
#endif
  }

  /**
   * @copydoc autopas::utils::Timer::stop()
   */
  long stop() {
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    return _timer.stop();
#else
    return 0;
#endif
  }

  /**
   * @copydoc autopas::utils::Timer::reset()
   */
  void reset() {
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    _timer.reset();
#endif
  }

  /**
   * @copydoc autopas::utils::Timer::addTime()
   */
  void addTime(long nanoseconds) {
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    _timer.addTime(nanoseconds);
#endif
  }

  /**
   * @copydoc autopas::utils::Timer::getTotalTime()
   */
  [[nodiscard]] long getTotalTime() const {
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    return _timer.getTotalTime();
#else
    return 0;
#endif
  }

 private:
  // The actual timer object is only instantiated if we are logging at TRACE level.
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
  Timer _timer;
#endif
};

}  // namespace autopas::utils