/**
 * @file EnergySensor.cpp
 * @author Maximilian Praus
 * @date 29.06.2024
 */

#include "autopas/utils/EnergySensor.h"

#include "autopas/options/EnergySensorOption.h"
#include "pmt.h"

namespace autopas::utils {

EnergySensor::EnergySensor(EnergySensorOption sensor) : _option(sensor) {
  if (_option != EnergySensorOption::none) {
    _sensor = pmt::Create(sensor.to_string());
  } else {
    AutoPasLog(WARN, "No sensor for energy consumption measurement specified. Energy will not be measured!");
  }
}

bool EnergySensor::startMeasurement() {
  if (_option != EnergySensorOption::none) {
    _start = _sensor->Read();
    return true;
  }
  return false;
}

const EnergySensorOption EnergySensor::getOption() const { return _option; }

bool EnergySensor::endMeasurement() {
  if (_option != EnergySensorOption::none) {
    _end = _sensor->Read();
    return true;
  }
  return false;
}

double EnergySensor::getJoules() const {
  if (_option != EnergySensorOption::none) {
    return _sensor->joules(_start, _end);
  }
  return -1;
}

double EnergySensor::getWatts() const {
  if (_option != EnergySensorOption::none) {
    return _sensor->watts(_start, _end);
  }
  return -1;
}

double EnergySensor::getSeconds() const {
  if (_option != EnergySensorOption::none) {
    return _sensor->seconds(_start, _end);
  }
  return -1;
}

long EnergySensor::getNanoJoules() const {
  if (_option != EnergySensorOption::none) {
    return static_cast<long>(_sensor->joules(_start, _end) * 1e9);
  }
  return -1;
}
}  // namespace autopas::utils