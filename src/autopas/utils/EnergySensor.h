/**
 * @file EnergySensor.h
 * @author Maximilian Praus
 * @date 29.06.2024
 */

#pragma once

#include "pmt.h"

namespace autopas::utils {


/**
 * Measure the energy consumption of a simulation on multiple hardwares.
 * Simulation needs to run on an Linux system to use Power measurement toolkit
 */

class EnergySensor {
    public:
    EnergySensor();

    bool startMeasurement();

    bool endMeasurement();

    double getJoules();

    double getWatts();

    double getSeconds();


private:
        std::unique_ptr<pmt::PMT> sensor; 
        pmt::State start;
        pmt::State end;

};

}