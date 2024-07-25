/**
 * @file EnergySensor.cpp
 * @author Maximilian Praus
 * @date 29.06.2024
 */

#include "autopas/utils/EnergySensor.h"
#include "autopas/options/EnergySensorOption.h"
#include "pmt.h"



namespace autopas::utils {

    EnergySensor::EnergySensor(EnergySensorOption sensor) : _sensor(pmt::Create(sensor.to_string())), _option(sensor) {
        if(_option == EnergySensorOption::none) {
        AutoPasLog(WARN, "No sensor for energy consumption measurement specified. Energy will not be measured!");
        }
    }

    bool EnergySensor::startMeasurement() {
        if(_option != EnergySensorOption::none) {
        _start = _sensor->Read();
        return true;
        } 
        return false;
    }

    const EnergySensorOption EnergySensor::getOption() {
        return _option;
    }

    bool EnergySensor::endMeasurement() {
        if(_option != EnergySensorOption::none) {
        _end = _sensor->Read();
        _energyTotal += _sensor->joules(_start, _end) * 1000;
        return true;
        }
        return false;
    }

    double EnergySensor::getJoules() {
        if (_option != EnergySensorOption::none) {
        return _sensor->joules(_start, _end);
        }
        return -1;
    }

    double EnergySensor::getWatts() {
        if(_option != EnergySensorOption::none) {
        return _sensor->watts(_start, _end);
        }
        return -1;
    }

    double EnergySensor::getSeconds() {
        if(_option != EnergySensorOption::none){
        return _sensor->seconds(_start, _end);
        }
        return -1;
    }

    long EnergySensor::getTotal() {
        return _energyTotal / 1000;
    }
}