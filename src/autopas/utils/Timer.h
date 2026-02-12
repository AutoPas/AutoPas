/**
 * @file Timer.h
 *
 * @date 11.01.2011
 * @author eckhardw
 */

#pragma once

#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>

namespace autopas::utils {

// Use steady_clock if high_resolution_clock is not monotonous in time.
using bestSteadyClock = std::conditional_t<std::chrono::high_resolution_clock::is_steady,
                                           std::chrono::high_resolution_clock, std::chrono::steady_clock>;

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
   * @return Total time in nano seconds.
   */
  [[nodiscard]] long getTotalTime() const { return _totalTime; }

  /**
   * Create a date stamp for the current moment with the given format.
   * @param format Date stamp format.
   * @return String representation of the current date in the given format.
   */
  static std::string getDateStamp(const std::string &format = "%Y-%m-%d_%H-%M-%S") {
    auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::ostringstream nowStrStr;
    tm unused;
    std::stringstream ss;
    ss << std::put_time(localtime_r(&now, &unused), format.c_str());
    return ss.str();
  }

 private:
  /**
   * Time point of last call of start().
   */
  bestSteadyClock::time_point _startTime;

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
