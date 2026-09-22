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
 * if AUTOPAS_ACTIVE_LEVEL is set to TRACE.
 * For a higher log level, functions are empty and will be removed by the compiler.
 */
class TraceTimer {
 public:
  /**
   * @copydoc autopas::utils::Timer::start()
   */
  void start() {
#if AUTOPAS_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    _timer.start();
#endif
  }

  /**
   * @copydoc autopas::utils::Timer::stop()
   */
  long stop() {
#if AUTOPAS_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    return _timer.stop();
#else
    return 0;
#endif
  }

  /**
   * @copydoc autopas::utils::Timer::reset()
   */
  void reset() {
#if AUTOPAS_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    _timer.reset();
#endif
  }

  /**
   * @copydoc autopas::utils::Timer::addTime()
   */
  void addTime(long nanoseconds) {
#if AUTOPAS_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    _timer.addTime(nanoseconds);
#endif
  }

  /**
   * @copydoc autopas::utils::Timer::getTotalTime()
   */
  [[nodiscard]] long getTotalTime() const {
#if AUTOPAS_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    return _timer.getTotalTime();
#else
    return 0;
#endif
  }

 private:
  // The actual timer object is only instantiated if we are logging at TRACE level.
#if AUTOPAS_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
  Timer _timer;
#endif
};

}  // namespace autopas::utils