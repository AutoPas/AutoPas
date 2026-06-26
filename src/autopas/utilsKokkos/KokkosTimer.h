/**
 * @file KokkosTimer.h
 * @date 17.06.2026
 * @author Luis Gall
 */

#pragma once

#ifdef AUTOPAS_ENABLE_KOKKOS

#include "Kokkos_Timer.hpp"

namespace autopas::utilsKokkos {

/**
 * Timer class to stop times.
 */
class KokkosTimer {
 public:
  KokkosTimer();

  virtual ~KokkosTimer();

  /**
   * start the timer.
   */
  void start();

  /**
   * Stops the timer and returns the time elapsed in nanoseconds since the last call to start.
   * It also adds the duration to the total time.
   * @return elapsed time in nanoseconds
   */
  long stop();

  /**
   * Resets the timer to 0.
   */
  void reset();

  /**
   * Adds the given amount of nanoseconds to the total time.
   * @param nanoseconds
   */
  void addTime(long nanoseconds);

  /**
   * Get total accumulated time.
   * @return Total time in nanoseconds.
   */
  [[nodiscard]] long getTotalTime() const { return _totalTime; }

 private:

  /**
   * Accumulated total time.
   */
  long _totalTime = 0;

  /**
   * Indicator if this timer currently is measuring.
   */
  bool _currentlyRunning = false;

  Kokkos::Timer _timer {};
};
}  // namespace autopas::utils

#endif