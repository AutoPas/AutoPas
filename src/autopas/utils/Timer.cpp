/**
 * @file Timer.cpp
 * @date 18.01.2011
 * @author tchipev
 */

#include "utils/Timer.h"
#include <iostream>

using namespace std::chrono;

autopas::utils::Timer::Timer() : _startTime{} {}

autopas::utils::Timer::~Timer() = default;

void autopas::utils::Timer::start() noexcept { _startTime = high_resolution_clock::now(); }

double autopas::utils::Timer::stop() {
  const high_resolution_clock::time_point time(high_resolution_clock::now());

  const duration<double> diff = duration_cast<duration<double>>(time - _startTime);

  return diff.count();
}
