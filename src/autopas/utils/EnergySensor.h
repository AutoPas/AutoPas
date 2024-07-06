/**
 * @file EnergySensor.h
 * @author Maximilian Praus
 * @date 29.06.2024
 */

#pragma once
#ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
#include "pmt.h"
#endif
#include "autopas/options/EnergySensorOption.h"

namespace autopas::utils {


/**
 * Measure the energy consumption of a simulation on multiple hardwares.
 * Simulation needs to run on an Linux system to use Power measurement toolkit
 */

class EnergySensor {
    public:
    EnergySensor(EnergySensorOption sensor);

    bool startMeasurement();

    bool endMeasurement();

    double getJoules();

    double getWatts();

    double getSeconds();

    /**
     * Getter for used sensor option
     * @return
     */
    const EnergySensorOption getOption();


private:
        #ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
        std::unique_ptr<pmt::PMT> _sensor; 
        pmt::State _start;
        pmt::State _end;
        #endif
        EnergySensorOption _option;

};

}