/*
 * Timer.cpp
 *
 * @Date: 11.01.2011
 * @Author: eckhardw
 */

#include "utils/Timer.h"
#include <iostream>

using namespace std;

utils::Timer::Timer() : _startTime{} {
  struct timespec info {};
  if (clock_getres(CLOCK_REALTIME, &info)) {
    std::cout << "Could not retrieve time resolution!" << endl;
  }
}

utils::Timer::~Timer() = default;

void utils::Timer::start() {
  if (clock_gettime(CLOCK_REALTIME, &_startTime)) {
    std::cout << "Could not retrieve time!" << endl;
  }
}

double utils::Timer::stop() {
  struct timespec time {};
  if (clock_gettime(CLOCK_REALTIME, &time)) {
    std::cout << "Could not retrieve time!" << endl;
  }

  double diff = time.tv_sec - _startTime.tv_sec;
  diff += ((double)(time.tv_nsec - _startTime.tv_nsec)) / 1000000000.0;

  return diff;
}
