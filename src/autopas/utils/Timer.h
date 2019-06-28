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
   * stops the timer and returns the time elapsed in seconds since the last call
   * to start
   * @return elapsed time in seconds
   */
  double stop();

 private:
  std::chrono::high_resolution_clock::time_point _startTime;
};
}  // namespace autopas
