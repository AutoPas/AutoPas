/**
 * @file EnergySensor.h
 * @author Maximilian Praus
 * @date 29.06.2024
 */

#pragma once
#include "autopas/options/EnergySensorOption.h"
#include "pmt.h"

namespace autopas::utils {

/**
 * Measure the energy consumption of a simulation on multiple hardwares.
 * Simulation needs to run on an Linux system to use Power measurement toolkit
 */

class EnergySensor {
 public:
  /**
   * Cosntructor for energy sensor. Takes @param sensor as possible energy sensors to use
   */
  EnergySensor(EnergySensorOption sensor);

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
  /**
   * Used energy sensor. Set in constuctor.
   */
  EnergySensorOption _option;
};

}  // namespace autopas::utils