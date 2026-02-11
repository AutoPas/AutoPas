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
   * Constructor for energy sensor. Takes @param sensor as possible energy sensors to use.
   * Causes a runtime error if the specified sensor is not available with pmt module used in AutoPas.
   * Currently, pmt in AutoPas can only be compiled with rapl and likwid.
   * Rapl sensor is the default sensor option. And it is currently forced to be compiled always.
   */
  EnergySensor(EnergySensorOption sensor);

  /**
   * Initializes the EnergySenor.
   * @param tuningMetricIsEnergy whether energy is selected as the tuning metric.
   * @return Bool whether energy measurement is turned on.
   */
  bool init(bool tuningMetricIsEnergy);

  /**
   * Start measurements. Sets _start
   * @return returns true if _sensor available. Otherwise false
   */
  bool startMeasurement();

  /**
   * Ending measurement. Sets _end
   * @return returns true if _sensor available. Otherwise false.
   */
  bool endMeasurement();

  /**
   * Get the energy consumed in Joules between the start and the end state
   * @return double for energy consumed in Joules.
   * Returns -1 when AutoPas in compiled without energy measurements enabled.
   */
  double getJoules() const;

  /**
   * Get the average power consumed in Watts between the start and the end state.
   * @return double for average power consumed in Watts.
   * Returns -1 when AutoPas in compiled without energy measurements enabled.
   */
  double getWatts() const;

  /**
   * Get seconds between current start and end state.
   * @return double for seconds between time stamps for energy measurements.
   * Returns -1 when AutoPas in compiled without energy measurements enabled.
   */
  double getEnergyDeltaT() const;

  /**
   * Getter for used sensor option.
   * @return EnergySensorOption displaying which energy sensor is used.
   * Returns -1 when AutoPas in compiled without energy measurements enabled.
   */
  const EnergySensorOption getOption() const;

  /**
   * Method to convert consumed energy in Joules to nanoJoules. This is used for tuning.
   * @return Energy consumed in nanoJoules.
   * Returns -1 when AutoPas in compiled without energy measurements enabled.
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