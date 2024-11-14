/**
 * @file EnergySensor.cpp
 * @author Maximilian Praus
 * @date 29.06.2024
 */

#include "autopas/utils/EnergySensor.h"

#include "autopas/options/EnergySensorOption.h"

namespace autopas::utils {

EnergySensor::EnergySensor(EnergySensorOption sensor) : _option(sensor) {
#ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
  _sensor = pmt::Create(sensor.to_string());
  // check if there is a way to check if sensor is created in PMT
#endif
}

bool EnergySensor::init(bool tuningMetricIsEnergy) {
#ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
  return true;
#else
  if (tuningMetricIsEnergy) {
    autopas::utils::ExceptionHandler::exception("Energy tuning cannot be performed with energy measurements disabled.");
  }
  return false;
#endif
}

bool EnergySensor::startMeasurement() {
#ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
  _start = _sensor->Read();
  return true;
#else
  return false;
#endif
}

const EnergySensorOption EnergySensor::getOption() const { return _option; }

bool EnergySensor::endMeasurement() {
#ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
  _end = _sensor->Read();
  return true;
#else
  return false;
#endif
}

double EnergySensor::getJoules() const {
#ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
  return _sensor->joules(_start, _end);
#else
  return -1;
#endif
}

double EnergySensor::getWatts() const {
#ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
  return _sensor->watts(_start, _end);
#else
  return -1;
#endif
}

double EnergySensor::getSeconds() const {
#ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
  return _sensor->seconds(_start, _end);
#else
  return -1;
#endif
}

long EnergySensor::getNanoJoules() const {
#ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
  return static_cast<long>(_sensor->joules(_start, _end) * 1e9);
#else
  return -1;
#endif
}
}  // namespace autopas::utils