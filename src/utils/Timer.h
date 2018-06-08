/**
 * @file Timer.h
 *
 * @date 11.01.2011
 * @author eckhardw
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <ctime>
namespace autopas {
namespace utils {

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
   * stops the timer and returns the time elapsed in seconds since the last call
   * to start
   * @return elapsed time in seconds
   */
  double stop();

 private:
  struct timespec _startTime;
};
}  // namespace utils
}  // namespace autopas
#endif /* TIMER_H_ */
