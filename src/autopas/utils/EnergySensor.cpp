/**
 * @file EnergySensor.cpp
 * @author Maximilian Praus
 * @date 29.06.2024
 */

#include "autopas/utils/EnergySensor.h"
#include "autopas/options/EnergySensorOption.h"

#ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
#include "pmt.h"
#endif



namespace autopas::utils {

    EnergySensor::EnergySensor(EnergySensorOption sensor) : _option(sensor) {
        #ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
        if(_option != EnergySensorOption::none) {
            AutoPasLog(WARN, "Sensor for energy measurement specified. Using {}", _option.to_string());
        this->_sensor = std::unique_ptr<pmt::PMT>(pmt::Create(sensor.to_string()));
        } else {
        AutoPasLog(WARN, "No sensor for energy consumption measurement specified. Energy will not be measured!");
        }
        #endif
    }

    bool EnergySensor::startMeasurement() {
        #ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
        if(_option != EnergySensorOption::none) {
        this->_start = this->_sensor->Read();
        return true;
        } 
        #endif
        return false;
    }

    const EnergySensorOption EnergySensor::getOption() {
        return _option;
    }

    bool EnergySensor::endMeasurement() {
        #ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
        if(_option != EnergySensorOption::none) {
        this->_end = this->_sensor->Read();
        return true;
        }
        #endif
        return false;
    }

    double EnergySensor::getJoules() {
        #ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
        if (_option != EnergySensorOption::none) {
        return this->_sensor->joules(this->_start, this->_end);
        }
        #endif
        return -1;
    }

    double EnergySensor::getWatts() {
        #ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
        if(_option != EnergySensorOption::none) {
        return this->_sensor->watts(this->_start, this->_end);
        }
        #endif
        return -1;
    }

    double EnergySensor::getSeconds() {
        #ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
        if(_option != EnergySensorOption::none){
        return this->_sensor->seconds(this->_start, this->_end);
        }
        #endif
        return -1;
    }
}