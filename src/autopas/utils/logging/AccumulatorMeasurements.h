/**
 * @file AccumulatorMeasurements.h
 * @author L. Llaveshi
 * @date 05.12.2024
 */

#pragma once

#include <string>

namespace autopas {

/**
 * Class to accumulate and store energy measurements across iterations.
 */
class AccumulatorMeasurements {
 public:
  /**
   * Constructor initializes all accumulated metrics to zero.
   */
  AccumulatorMeasurements();

  /**
   * Adds energy measurements to the accumulator.
   * @param psys Energy consumed by the entire system (Psys).
   * @param pkg Energy consumed by the CPU package (Pkg).
   * @param ram Energy consumed by the RAM.
   */
  void addEnergy(double psys, double pkg, double ram);

  /**
   * Gets the total accumulated system energy.
   * @return Total accumulated energy (Psys).
   */
  double getAccumulatedEnergyPsys() const;

  /**
   * Gets the total accumulated package energy.
   * @return Total accumulated energy (Pkg).
   */
  double getAccumulatedEnergyPkg() const;

  /**
   * Gets the total accumulated RAM energy.
   * @return Total accumulated energy (Ram).
   */
  double getAccumulatedEnergyRam() const;

  /**
   * Gets a formatted string representation of all accumulated measurements.
   * @return Formatted string.
   */
  std::string toString() const;

 private:
  double accumulatedEnergyPsys;   ///< Accumulated system energy (Psys).
  double accumulatedEnergyPkg;    ///< Accumulated package energy (Pkg).
  double accumulatedEnergyRam;    ///< Accumulated RAM energy.
};

}  // namespace autopas
