/**
 * @file Timer.cpp
 * @date 18.01.2011
 * @author tchipev
 */

#include "utils/Timer.h"

#include "ExceptionHandler.h"

using namespace std::chrono;

autopas::utils::Timer::Timer() : _startTime{} {}

autopas::utils::Timer::~Timer() = default;

void autopas::utils::Timer::start() {
  if (_currentlyRunning) {
    ExceptionHandler::exception("Trying to start a timer that is already started!");
  }
  _currentlyRunning = true;
  _startTime = high_resolution_clock::now();
}

long autopas::utils::Timer::stop() {
  const auto time(high_resolution_clock::now());

  if (not _currentlyRunning) {
    ExceptionHandler::exception("Trying to stop a timer that was not started!");
  }
  _currentlyRunning = false;

  const auto diff = duration_cast<nanoseconds>(time - _startTime).count();

  _totalTime += diff;

  return _lapTime;
}

void autopas::utils::Timer::reset() {
  _totalTime = 0;
  _lapTime = 0;
}

void autopas::utils::Timer::resetLap() { _lapTime = 0; }

void autopas::utils::Timer::addTime(long nanoseconds) {
  _totalTime += nanoseconds;
  _lapTime += nanoseconds;
}
