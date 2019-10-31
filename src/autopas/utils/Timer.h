/**
 * @file Timer.h
 *
 * @date 11.01.2011
 * @author eckhardw
 */

#pragma once

#include <chrono>

namespace autopas::utils {

/**
 * Timer class to stop times.
 */
class Timer {
 public:
  Timer();

  virtual ~Timer();

  /**
   * start the timer.
   */
  void start() noexcept;

  /**
   * Stops the timer and returns the time elapsed in seconds since the last call
   * to start and adds the duration to the total time.
   * @return elapsed time in seconds
   */
  double stop();

  /**
   * Get total accumulated time.
   * @return
   */
  long getTotalTime() const { return _totalTime; }

 private:
  /**
   * Time point of last call of start().
   */
  std::chrono::high_resolution_clock::time_point _startTime;

  /**
   * Accumulated total time.
   */
  long _totalTime = 0;

  /**
   * Indicator if this timer currently is measuring.
   */
   bool _currentlyRunning = false;
};
}  // namespace autopas::utils
