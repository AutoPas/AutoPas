/**
 * @file EnergySensor.cpp
 * @author Maximilian Praus
 * @date 29.06.2024
 */

#include "autopas/utils/EnergySensor.h"

#include "pmt.h"



namespace autopas::utils {

    EnergySensor::EnergySensor() {
        this->sensor = std::unique_ptr<pmt::PMT>(pmt::Create("rapl"));
    }

    bool EnergySensor::startMeasurement() {
        this->start = this->sensor->Read();
        return true;
    }

    bool EnergySensor::endMeasurement() {
        this->end = this->sensor->Read();
        return true;
    }

    double EnergySensor::getJoules() {
        return this->sensor->joules(this->start, this->end);
    }

    double EnergySensor::getWatts() {
        return this->sensor->watts(this->start, this->end);
    }

    double EnergySensor::getSeconds() {
        return this->sensor->seconds(this->start, this->end);
    }
}