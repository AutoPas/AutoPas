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

/**
 * Timer class to measure timings.
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
   * It also adds the duration to the total time and the lap time.
   * @return elapsed time in nanoseconds
   */
  long stop();

  /**
   * Resets the timer (and lap timer) to 0.
   */
  void reset();

  /**
   * Resets "lap" counters to 0.
   */
  void resetLap();

  /**
   * Adds the given amount of nanoseconds to the total time and lap time.
   * @param nanoseconds
   */
  void addTime(long nanoseconds);

  /**
   * Get total accumulated time.
   * @return Total time in nanoseconds.
   */
  [[nodiscard]] long getTotalTime() const { return _totalTime; }

  /**
   * Get the total accumulated time since the start of the "lap".
   * @return Time in nanoseconds.
   */
  long getLapTime() const { return _lapTime; }

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
  std::chrono::high_resolution_clock::time_point _startTime;

  /**
   * Accumulated total time.
   */
  long _totalTime = 0;

  /**
   * Accumulated "lap" time. This can serve as an additional timer accumulator which can be reset without affecting
   * the total time, similarly to laps in many stopwatches.
   */
  long _lapTime = 0;

  /**
   * Indicator if this timer currently is measuring.
   */
  bool _currentlyRunning = false;
};
}  // namespace autopas::utils
