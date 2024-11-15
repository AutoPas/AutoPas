/**
 * @file EnergySensor.h
 * @author Maximilian Praus
 * @date 29.06.2024
 */

#pragma once
#include "autopas/options/EnergySensorOption.h"
#ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
#include "pmt.h"
#endif

namespace autopas::utils {

/**
 * Measure the energy consumption of a simulation on multiple hardwares.
 * Simulation needs to run on an Linux system to use Power measurement toolkit
 */

class EnergySensor {
 public:
  /**
   * Constructor for energy sensor. Takes @param sensor as possible energy sensors to use
   */
  EnergySensor(EnergySensorOption sensor);

  /**
   * Initializes the EnergySenor
   * @param tuningMetricIsEnergy whether energy is selected as the tuning metric
   * @return whether energy measurement is turned on
   */
  bool init(bool tuningMetricIsEnergy);

  /**
   * Start measurements. Sets _start
   * @return returns true if _sensor available. Otherwise false
   */
  bool startMeasurement();

  /**
   * Ending measurement. Sets _end
   * @return returns true if _sensor available. Otherwise false
   */
  bool endMeasurement();

  /**
   * Get joules consumed between start and end state
   * @return double for watts consumed
   */
  double getJoules() const;

  /**
   * Get watts consumed between start and end state
   * @return double for watts consumed
   */
  double getWatts() const;

  /**
   * Get seconds between current start and end state
   * @return double for seconds between time stamps
   */
  double getSeconds() const;

  /**
   * Getter for used sensor option
   * @return EnergySensorOption displaying which energy sensor is used
   */
  const EnergySensorOption getOption() const;

  /**
   * Method to convert consumed joules to nanojoules. Used for tuning
   * @return consumed nanojoules
   */
  long getNanoJoules() const;

 private:
#ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
  /**
   * Pointer to the pmt instance to read the energy from
   */
  std::unique_ptr<pmt::PMT> _sensor;
  /**
   * Current start state. Marks beginning of current energy cycle
   */
  pmt::State _start;
  /**
   * Current end state. Marks ending of current energy measurement cycle
   */
  pmt::State _end;
#endif
  /**
   * Used energy sensor. Set in constuctor.
   */
  EnergySensorOption _option;
};

}  // namespace autopas::utils