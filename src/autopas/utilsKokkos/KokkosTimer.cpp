/**
 * @file KokkosTimer.cpp
 * @date 17.06.2026
 * @author Luis Gall
 */

#ifdef AUTOPAS_ENABLE_KOKKOS

#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utilsKokkos/KokkosTimer.h"

autopas::utilsKokkos::KokkosTimer::KokkosTimer() {}

autopas::utilsKokkos::KokkosTimer::~KokkosTimer() = default;

void autopas::utilsKokkos::KokkosTimer::start() {
  if (_currentlyRunning) {
    utils::ExceptionHandler::exception("Trying to start a timer that is already started!");
  }
  _currentlyRunning = true;

  _timer.reset();
}

long autopas::utilsKokkos::KokkosTimer::stop() {

  if (not _currentlyRunning) {
    utils::ExceptionHandler::exception("Trying to stop a timer that was not started!");
  }
  _currentlyRunning = false;

  auto elapsedTime = _timer.seconds() * 1e9;
  _totalTime += static_cast<long>(elapsedTime);

  return elapsedTime;
}

void autopas::utilsKokkos::KokkosTimer::reset() { _totalTime = 0; }

void autopas::utilsKokkos::KokkosTimer::addTime(long nanoseconds) { _totalTime += nanoseconds; }

#endif